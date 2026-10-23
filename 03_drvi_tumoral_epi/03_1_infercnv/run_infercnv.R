#!/usr/bin/env Rscript

# 03_1 run_infercnv: inferCNV on one cohort, reduced to a per-cell table.
#
# Reference = ref_tcell + ref_myeloid; subclusters are computed by leiden inside the run.
# inferCNV's working directory (1-3 GB per cohort) is deleted on success. The malignant
# call is made afterwards, in call_malignant.ipynb.
#
# Input : $DATA_DIR/03_tum/input/<cohort>/        (prepare_infercnv_input.py)
#         $DATA_DIR/03_tum/gene_order_hg38_gencode_v27.txt
# Output: $DATA_DIR/03_tum/summary/<cohort>_cnv.csv, one row per cell:
#           group        ref_tcell / ref_myeloid / stromal / immune_heldout / epi
#           subcluster   per-run leiden subcluster "<group>_s<i>", the unit of the verdict
#           cnv_score    mean squared residual against the mean of the reference cells
#           cnv_corr     correlation with the mean residual of the top 5% epithelial cells
#           chr1..chr22  mean residual per chromosome
#         $FIG_DIR/infercnv_<cohort>.png      (default ../figures/03_1_infercnv/)
#
# Usage:
#   conda activate infercnv-r
#   export DATA_DIR=<data root>
#   Rscript run_infercnv.R --cohort Patient64
#   Rscript run_infercnv.R --cohort Patient64 --threads 8 --keep-work
#   Rscript run_infercnv.R --cohort Patient64 --hmm    # also the i6 HMM (hours, needs JAGS)
#
# infercnv_all.sh runs it on every prepared cohort.

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(infercnv)
})

# Needed with subclusters: inferCNV passes the trees through Newick and cannot parse
# branch lengths printed as "1e-05".
options(scipen = 100)

option_list <- list(
  make_option("--cohort", type = "character", help = "cohort to run (required)"),
  make_option("--data-dir", type = "character", default = Sys.getenv("DATA_DIR"),
              dest = "data_dir", help = "datasets directory [env DATA_DIR]"),
  make_option("--fig-dir", type = "character", default = Sys.getenv("FIG_DIR"),
              dest = "fig_dir", help = "figure directory [env FIG_DIR]"),
  make_option("--threads", type = "integer", default = 8L,
              help = "threads for infercnv::run() [%default]"),
  make_option("--hmm", action = "store_true", default = FALSE,
              help = "also run the i6 HMM (hours per cohort; not needed for the binary call)"),
  make_option("--analysis-mode", type = "character", default = "subclusters",
              dest = "analysis_mode",
              help = "'subclusters' (default) or 'samples' to reproduce the pre-subcluster runs"),
  make_option("--leiden-resolution", type = "character", default = "0.005",
              dest = "leiden_resolution",
              help = "leiden resolution for the subclustering, or 'auto' [%default]"),
  make_option("--keep-work", action = "store_true", default = FALSE, dest = "keep_work",
              help = "keep inferCNV's working directory (1-3 GB per cohort)"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "re-run a cohort whose summary already exists")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cohort)) stop("--cohort is required")
if (!nzchar(opt$data_dir)) stop("set DATA_DIR (or pass --data-dir)")

# Parameters
REF_GROUPS <- c("ref_tcell", "ref_myeloid")
CUTOFF <- 0.1 # inferCNV's recommended value for 10x data
TOP_FRAC <- 0.05 # top epithelial cells defining the cohort's CNV profile for cnv_corr
MIN_TOP_CELLS <- 20 # but never fewer than this
GENE_ORDER_NAME <- "gene_order_hg38_gencode_v27.txt"
set.seed(0)

cohort   <- opt$cohort
data_dir <- normalizePath(path.expand(opt$data_dir), mustWork = TRUE)
cnv_dir  <- file.path(data_dir, "03_tum")
in_dir   <- file.path(cnv_dir, "input", cohort)
work_dir <- file.path(cnv_dir, "work", cohort)
sum_dir  <- file.path(cnv_dir, "summary")
sum_path <- file.path(sum_dir, paste0(cohort, "_cnv.csv"))

fig_dir <- if (nzchar(opt$fig_dir)) opt$fig_dir else
  file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))),
            "figures", "03_1_infercnv")
# infercnv::run() does not create parent directories.
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sum_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (file.exists(sum_path) && !opt$force) {
  cat(sprintf("[have] %s: %s already exists, skipping (--force to re-run)\n", cohort, sum_path))
  quit(status = 0)
}
if (!dir.exists(in_dir)) stop(sprintf("no prepared input at %s; run prepare_infercnv_input.py", in_dir))

cat(sprintf(paste0("cohort   : %s\ninput    : %s\nwork dir : %s\nfig dir  : %s\n",
                   "threads  : %d | HMM: %s | mode: %s | leiden res: %s\n\n"),
            cohort, in_dir, work_dir, fig_dir, opt$threads, opt$hmm,
            opt$analysis_mode, opt$leiden_resolution))

# Load (Matrix Market: tens of MB per cohort instead of ~1 GB as a dense table)
counts <- readMM(file.path(in_dir, "counts.mtx"))
counts <- as(as(counts, "CsparseMatrix"), "dgCMatrix")
rownames(counts) <- readLines(file.path(in_dir, "genes.tsv"))
colnames(counts) <- readLines(file.path(in_dir, "barcodes.tsv"))
cat(sprintf("counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))

annot <- read.delim(file.path(in_dir, "annotations.tsv"), header = FALSE,
                    col.names = c("cell", "group"), stringsAsFactors = FALSE)
print(table(annot$group))
stopifnot(all(REF_GROUPS %in% annot$group))

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = file.path(in_dir, "annotations.tsv"),
  gene_order_file   = file.path(cnv_dir, GENE_ORDER_NAME),
  ref_group_names   = REF_GROUPS
)
rm(counts); invisible(gc())

# Run
#   denoise            the summary uses the denoised residuals
#   cluster_by_groups  heatmap shows each group separately
#   analysis_mode      "subclusters": per-patient leiden partition, the unit of the verdict
#   leiden_resolution  fixed, because "auto" shrinks as the cohort grows
#   no_prelim_plot, write_expr_matrix, save_rds, save_final_rds, resume_mode: skip GB-scale files
t0 <- Sys.time()
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff            = CUTOFF,
  out_dir           = work_dir,
  cluster_by_groups = TRUE,
  denoise           = TRUE,
  HMM               = opt$hmm,
  HMM_type          = "i6",
  analysis_mode     = opt$analysis_mode,
  tumor_subcluster_partition_method = "leiden",
  leiden_resolution = if (opt$leiden_resolution == "auto") "auto" else
                        as.numeric(opt$leiden_resolution),
  num_threads       = opt$threads,
  no_prelim_plot    = TRUE,
  output_format     = "png",
  write_expr_matrix = FALSE,
  save_rds          = FALSE,
  save_final_rds    = FALSE,
  resume_mode       = FALSE
)
cat(sprintf("\ninfercnv::run() took %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Per-cell table
expr <- infercnv_obj@expr.data
ref_idx <- sort(unlist(infercnv_obj@reference_grouped_cell_indices, use.names = FALSE))
stopifnot(length(ref_idx) > 0)

# Residual against the reference mean, so inferCNV's own centring does not enter the score.
ref_mean <- rowMeans(expr[, ref_idx, drop = FALSE])
resid <- expr - ref_mean
rm(expr); invisible(gc())

cnv_score <- colMeans(resid^2)

group <- setNames(annot$group, annot$cell)[colnames(resid)]
stopifnot(!anyNA(group))

# Subclusters live only in the returned object: "<group>_s<i>" -> column indices
# (cell names in older inferCNV versions).
subcluster <- rep(NA_character_, ncol(resid))
sc <- infercnv_obj@tumor_subclusters$subclusters
if (is.null(sc)) {
  cat("no subclusters in the object (analysis_mode='samples'?)\n")
} else {
  for (nm in names(sc)) {
    v <- sc[[nm]]
    if (is.list(v)) { # nested one level in some versions
      for (nm2 in names(v)) {
        ix <- v[[nm2]]
        if (is.character(ix)) ix <- match(ix, colnames(resid))
        subcluster[ix] <- nm2
      }
    } else {
      if (is.character(v)) v <- match(v, colnames(resid))
      subcluster[v] <- nm
    }
  }
  cat(sprintf("subclusters: %d over %d cells (%d cells unassigned)\n",
              length(unique(na.omit(subcluster))), ncol(resid), sum(is.na(subcluster))))
  print(table(group, ifelse(is.na(subcluster), "<none>", "assigned")))
}

# cnv_corr: correlation with the mean profile of the most aneuploid epithelial cells.
epi_cells <- which(group == "epi")
n_top <- max(MIN_TOP_CELLS, ceiling(TOP_FRAC * length(epi_cells)))
n_top <- min(n_top, length(epi_cells))
top_cells <- epi_cells[order(cnv_score[epi_cells], decreasing = TRUE)[seq_len(n_top)]]
top_profile <- rowMeans(resid[, top_cells, drop = FALSE])
cnv_corr <- as.vector(cor(resid, top_profile))
cat(sprintf("CNV reference profile built on the top %d of %d epithelial cells\n",
            n_top, length(epi_cells)))

# Mean residual per chromosome, to see which chromosomes drive a call.
chrs <- infercnv_obj@gene_order$chr[match(rownames(resid), rownames(infercnv_obj@gene_order))]
chr_levels <- paste0("chr", 1:22)
chr_means <- sapply(chr_levels, function(ch) {
  rows <- which(as.character(chrs) == ch)
  if (length(rows) == 0) rep(NA_real_, ncol(resid)) else colMeans(resid[rows, , drop = FALSE])
})

out <- data.frame(
  cell       = colnames(resid),
  cohort     = cohort,
  group      = unname(group),
  subcluster = subcluster,
  cnv_score  = unname(cnv_score),
  cnv_corr   = cnv_corr,
  chr_means,
  check.names = FALSE, stringsAsFactors = FALSE
)
write.csv(out, sum_path, row.names = FALSE)
cat(sprintf("Wrote %s (%d cells)\n", sum_path, nrow(out)))

cat("\nmedian cnv_score by group:\n")
print(round(tapply(out$cnv_score, out$group, median), 5))

# Quick look: epithelial subclusters against the stromal median.
epi_sub <- out[out$group == "epi" & !is.na(out$subcluster), ]
if (nrow(epi_sub)) {
  cat("\nepithelial subclusters (n cells, median cnv_score):\n")
  print(round(do.call(rbind, lapply(split(epi_sub, epi_sub$subcluster), function(d)
    data.frame(n = nrow(d), median_cnv_score = median(d$cnv_score)))), 5))
  cat(sprintf("stromal median for reference: %.5f\n",
              median(out$cnv_score[out$group == "stromal"])))
}

# Keep the heatmap, delete the working directory
heatmap_src <- file.path(work_dir, "infercnv.png")
if (file.exists(heatmap_src)) {
  file.copy(heatmap_src, file.path(fig_dir, sprintf("infercnv_%s.png", cohort)), overwrite = TRUE)
  cat(sprintf("Wrote %s\n", file.path(fig_dir, sprintf("infercnv_%s.png", cohort))))
} else {
  warning(sprintf("no heatmap at %s", heatmap_src))
}

if (opt$keep_work) {
  cat(sprintf("[keep] working directory left at %s\n", work_dir))
} else {
  unlink(work_dir, recursive = TRUE)
  cat(sprintf("[clean] removed %s\n", work_dir))
}
