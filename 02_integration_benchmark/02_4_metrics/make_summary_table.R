#!/usr/bin/env Rscript
#
# 02_4 metrics: final summary table, methods ranked by
#     overall = 0.6 * bio conservation + 0.4 * batch correction
# on min-max scaled metrics (Luecken et al. 2022). Draws the published scIB figure
# with the vendored code in ../utils/, or falls back to a plain heatmap.
#
# Usage:
#   Rscript make_summary_table.R -i <merged.csv>
#   Rscript make_summary_table.R -i <merged.csv> -o <outdir> --viz-dir <path>
#
# After the metrics have run on the cluster:
#   export DATA_DIR=<data root>
#   rsync -av $CLUSTER_HOST:$CLUSTER_DATA_DIR/02_metrics/ $DATA_DIR/02_metrics/
#   cd 02_integration_benchmark/02_4_metrics
#   conda activate benchmark-py-r
#   python3 merge_metrics.py -o $DATA_DIR/02_metrics_merged.csv -r $DATA_DIR/02_metrics \
#       --glob "$DATA_DIR/02_metrics/shiao/metrics/*/hvg/*.csv"
#   Rscript make_summary_table.R -i $DATA_DIR/02_metrics_merged.csv
# (default output: the phase's figures/02_4_metrics/)

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "merged metrics CSV"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "output directory [default: the phase's figures/02_4_metrics]"),
  make_option("--viz-dir", type = "character", default = NULL,
              dest = "viz_dir", help = "directory with plotSingleTaskRNA.R + knit_table.R + img/"),
  make_option("--weight-batch", type = "double", default = 0.4, dest = "weight_batch",
              help = "weight of the batch-correction score [default %default]")
)))
stopifnot(!is.null(opt$input))

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") script_dir <- "."
viz_dir <- if (!is.null(opt$viz_dir)) opt$viz_dir else file.path(script_dir, "..", "utils")

# Default anchored to the script, not the working directory.
out_dir <- if (!is.null(opt$output)) opt$output else
  file.path(script_dir, "..", "figures", "02_4_metrics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
opt$output <- out_dir

# Metric groups for the fallback scorer.
BATCH_METRICS <- c("PCR_batch", "ASW_label/batch", "iLISI", "graph_conn")
BIO_METRICS   <- c("NMI_cluster/label", "ARI_cluster/label", "ASW_label",
                   "isolated_label_F1", "isolated_label_silhouette", "cLISI",
                   "cell_cycle_conservation", "hvg_overlap")


# the published figure, via the vendored scIB code 
run_scib_plot <- function() {
  plotter <- file.path(viz_dir, "plotSingleTaskRNA.R")
  knit    <- file.path(viz_dir, "knit_table.R")
  if (!file.exists(plotter) || !file.exists(knit)) return(FALSE)

  # Packages plotSingleTaskRNA.R needs at source time.
  needed <- c("tibble", "RColorBrewer", "dynutils", "stringr", "Hmisc", "plyr",
              "ggplot2", "cowplot", "dplyr", "ggimage", "scales", "ggforce")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("[summary] vendored scIB plotting needs these R packages: ",
            paste(missing, collapse = ", "),
            "\n           install with: Rscript -e \"install.packages(c(",
            paste(sprintf("'%s'", missing), collapse = ", "), "))\"")
    return(FALSE)
  }

  # The scIB parser splits <method>_<type> on "_", so strip underscores from
  # method names (seurat_cca, seurat_rpca) in the plotter's copy only.
  fixed <- read.csv(opt$input, check.names = FALSE)
  fixed[[1]] <- vapply(as.character(fixed[[1]]), function(lbl) {
    parts <- strsplit(lbl, "/")[[1]]
    toks <- strsplit(parts[length(parts)], "_")[[1]]
    type <- toks[length(toks)]
    method <- paste(toks[-length(toks)], collapse = "_")
    method <- switch(method,
                     "seurat_cca" = "seurat",       # -> "Seurat v3 CCA"
                     "seurat_rpca" = "seuratrpca",  # -> "Seurat v3 RPCA"
                     gsub("_", "", method))
    parts[length(parts)] <- paste0(method, "_", type)
    paste(parts, collapse = "/")
  }, character(1))

  in_abs <- file.path(normalizePath(opt$output), ".plot_input.csv")
  write.csv(fixed, in_abs, row.names = FALSE)

  out_abs <- normalizePath(opt$output)
  old_wd  <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  # The vendored code reads knit_table.R and img/ relative to the working dir.
  setwd(viz_dir)
  source("plotSingleTaskRNA.R")
  plotSingleTaskRNA(csv_metrics_path = in_abs, outdir = out_abs,
                    weight_batch = opt$weight_batch)
  message("[summary] wrote the published-style figure(s) to ", out_abs)
  TRUE
}

# fallback: self-contained min-max scorer + plain heatmap 
run_fallback <- function() {
  suppressPackageStartupMessages(library(ggplot2))
  message("[summary] using the built-in scorer (plain heatmap, not the published look).")

  tab <- read.csv(opt$input, row.names = 1, check.names = FALSE)
  run_label <- sub(".*/", "", rownames(tab))

  present <- function(cols) cols[cols %in% colnames(tab)]
  batch_cols <- present(BATCH_METRICS)
  bio_cols   <- present(BIO_METRICS)

  minmax <- function(x) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1]) || diff(rng) == 0) return(rep(NA_real_, length(x)))
    (x - rng[1]) / diff(rng)
  }
  scaled <- as.data.frame(lapply(tab, minmax), check.names = FALSE)
  rownames(scaled) <- rownames(tab)

  row_mean <- function(df, cols) {
    if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
    rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  bio_score   <- row_mean(scaled, bio_cols)
  batch_score <- row_mean(scaled, batch_cols)
  overall <- (1 - opt$weight_batch) * bio_score + opt$weight_batch * batch_score

  scores <- data.frame(run = run_label,
                       bio_conservation = round(bio_score, 4),
                       batch_correction = round(batch_score, 4),
                       overall = round(overall, 4), check.names = FALSE)
  scores <- scores[order(-scores$overall), ]
  write.csv(scores, file.path(opt$output, "summary_scores.csv"), row.names = FALSE)

  long <- data.frame(
    run = factor(rep(run_label, ncol(scaled)), levels = rev(scores$run)),
    metric = factor(rep(colnames(scaled), each = nrow(scaled)), levels = c(bio_cols, batch_cols)),
    value = as.vector(as.matrix(scaled)))
  long <- long[long$metric %in% c(bio_cols, batch_cols), ]

  p <- ggplot(long, aes(metric, run, fill = value)) +
    geom_tile(colour = "grey90") +
    scale_fill_viridis_c(na.value = "grey95", limits = c(0, 1)) +
    labs(x = NULL, y = NULL, fill = "min-max\nscaled",
         title = "Integration benchmark (overall = 0.6 bio + 0.4 batch)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(opt$output, "summary_heatmap.png"), p,
         width = 9, height = 0.4 * nrow(scores) + 2, dpi = 150)
  message("[summary] wrote summary_scores.csv + summary_heatmap.png to ", opt$output)
}


if (!run_scib_plot()) run_fallback()
