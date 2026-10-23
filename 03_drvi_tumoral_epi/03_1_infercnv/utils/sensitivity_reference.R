#!/usr/bin/env Rscript

# 03_1 sensitivity_reference: inferCNV on one reference configuration.
#
# Same call and parameters as run_infercnv.R, but the reference groups are read from the
# ref_groups.txt written by sensitivity_reference.py --mode prepare.
#
# Input : <config-dir>/{counts.mtx, genes.tsv, barcodes.tsv, annotations.tsv, ref_groups.txt}
# Output: <config-dir>_cells.csv (cell, group, cnv_score, cnv_corr); working directory deleted
#
# Usage:
#   conda activate infercnv-r
#   export DATA_DIR=<data root>
#   Rscript utils/sensitivity_reference.R --config-dir $DATA_DIR/03_tum/sensitivity/reference/Patient16/A_current

suppressPackageStartupMessages({library(optparse); library(Matrix); library(infercnv)})
options(scipen = 100) # see run_infercnv.R
set.seed(0)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config-dir", type = "character", dest = "config_dir"),
  make_option("--threads", type = "integer", default = 20L),
  make_option("--resolution", type = "character", default = "0.005")
)))
stopifnot(!is.null(opt$config_dir))

d <- normalizePath(opt$config_dir, mustWork = TRUE)
data_dir <- dirname(dirname(dirname(dirname(d)))) # .../03_tum/sensitivity/reference/<coh>/<cfg>
gene_order <- file.path(dirname(dirname(dirname(d))), "..", "gene_order_hg38_gencode_v27.txt")
gene_order <- normalizePath(gene_order, mustWork = TRUE)
out_csv <- file.path(dirname(d), paste0(basename(d), "_cells.csv"))
if (file.exists(out_csv)) { cat(sprintf("[have] %s\n", out_csv)); quit(status = 0) }

counts <- as(as(readMM(file.path(d, "counts.mtx")), "CsparseMatrix"), "dgCMatrix")
rownames(counts) <- readLines(file.path(d, "genes.tsv"))
colnames(counts) <- readLines(file.path(d, "barcodes.tsv"))
refs <- readLines(file.path(d, "ref_groups.txt"))
annot <- read.delim(file.path(d, "annotations.tsv"), header = FALSE,
                    col.names = c("cell", "group"), stringsAsFactors = FALSE)
cat(sprintf("%s / %s: %d genes x %d cells | reference: %s\n", basename(dirname(d)),
            basename(d), nrow(counts), ncol(counts), paste(refs, collapse = ", ")))

work <- file.path(d, "work")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
obj <- CreateInfercnvObject(raw_counts_matrix = counts,
                            annotations_file = file.path(d, "annotations.tsv"),
                            gene_order_file = gene_order, ref_group_names = refs)
rm(counts); invisible(gc())
obj <- infercnv::run(obj, cutoff = 0.1, out_dir = work, cluster_by_groups = TRUE,
                     denoise = TRUE, HMM = FALSE, analysis_mode = "subclusters",
                     tumor_subcluster_partition_method = "leiden",
                     leiden_resolution = if (opt$resolution == "auto") "auto"
                                         else as.numeric(opt$resolution),
                     num_threads = opt$threads, no_prelim_plot = TRUE, output_format = "png",
                     write_expr_matrix = FALSE, save_rds = FALSE, save_final_rds = FALSE,
                     resume_mode = FALSE)

expr <- obj@expr.data
ridx <- sort(unlist(obj@reference_grouped_cell_indices, use.names = FALSE))
resid <- expr - rowMeans(expr[, ridx, drop = FALSE])
rm(expr); invisible(gc())
score <- colMeans(resid^2)
group <- setNames(annot$group, annot$cell)[colnames(resid)]
ei <- which(group == "epi")
nt <- min(length(ei), max(20, ceiling(0.05 * length(ei))))
tp <- rowMeans(resid[, ei[order(score[ei], decreasing = TRUE)[seq_len(nt)]], drop = FALSE])

write.csv(data.frame(cell = colnames(resid), group = unname(group),
                     cnv_score = unname(score), cnv_corr = as.vector(cor(resid, tp)),
                     stringsAsFactors = FALSE), out_csv, row.names = FALSE)
cat(sprintf("Wrote %s\n", out_csv))
unlink(work, recursive = TRUE)
