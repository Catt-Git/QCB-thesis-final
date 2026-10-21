#!/usr/bin/env Rscript

# 02_2 integration: dispatcher for the R methods (fastmnn, seurat_cca, seurat_rpca;
# bodies in integration_methods.R). Writes the integrated Seurat .rds; run_all.sh
# then converts it with rds_to_h5ad.py. The HVG vector is only used as Seurat
# anchor features (the input is already HVG-only).
#
# Usage:
#   export DATA_DIR=<data root>
#   Rscript run_integration.R -m seurat_cca \
#       -i $DATA_DIR/shiao_hvg_2k.rds \
#       -o $DATA_DIR/02_integration/seurat_cca_unscaled.rds \
#       -v $DATA_DIR/shiao_hvg_2k_unintegrated_list.rds

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
})

# Source the method bodies from this script's directory.
.args_all <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("^--file=", "", .args_all[grep("^--file=", .args_all)])
.script_dir <- if (length(.script_path)) dirname(normalizePath(.script_path)) else getwd()
source(file.path(.script_dir, "integration_methods.R"))

option_list <- list(
  make_option(c("-m", "--method"), type = "character",
              help = "one of: fastmnn, seurat_cca, seurat_rpca"),
  make_option(c("-i", "--input"), type = "character",
              help = "input Seurat .rds (from 02_1 h5ad_to_rds.R)"),
  make_option(c("-o", "--output"), type = "character",
              help = "output integrated Seurat .rds"),
  make_option(c("-b", "--batch"), type = "character", default = "cohort",
              help = "batch column in meta.data [default %default]"),
  make_option(c("-v", "--hvg"), type = "character", default = NA,
              help = "HVG list .rds (required for the Seurat methods)"),
  make_option("--reference", type = "character",
              default = "Patient53,Patient16,Patient43",
              help = "comma-separated reference batches for CCA/RPCA [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$method) || is.null(opt$input) || is.null(opt$output)) {
  stop("-m/--method, -i/--input and -o/--output are required", call. = FALSE)
}
if (!opt$method %in% names(INTEGRATION_METHODS)) {
  stop(sprintf("unknown method %s; available: %s", opt$method,
               paste(names(INTEGRATION_METHODS), collapse = ", ")), call. = FALSE)
}
if (!file.exists(opt$input)) {
  stop(sprintf("input not found: %s", opt$input), call. = FALSE)
}

reference_batches <- trimws(strsplit(opt$reference, ",")[[1]])

# The Seurat methods need the HVG vector; fastMNN does not.
hvg <- NULL
if (!is.na(opt$hvg)) {
  hvg <- unlist(readRDS(opt$hvg), use.names = FALSE)
} else if (opt$method %in% c("seurat_cca", "seurat_rpca")) {
  stop(sprintf("method %s needs -v/--hvg (the HVG list .rds)", opt$method),
       call. = FALSE)
}

message(sprintf("[read] %s", opt$input))
sobj <- readRDS(opt$input)
cells_before <- colnames(sobj)
message(sprintf("[read] %d genes x %d cells, %d batches in '%s'",
                nrow(sobj), ncol(sobj),
                length(unique(sobj@meta.data[[opt$batch]])), opt$batch))

message(sprintf("[run] %s", opt$method))
integrated <- INTEGRATION_METHODS[[opt$method]](sobj, opt$batch, hvg, reference_batches)

# Metrics need the same cells in the same order; IntegrateData can drop some.
if (ncol(integrated) != length(cells_before)) {
  stop(sprintf("cell count changed during integration: %d -> %d",
               length(cells_before), ncol(integrated)), call. = FALSE)
}
if (!identical(colnames(integrated), cells_before)) {
  # Reorder to the input order if the cells are the same set but permuted.
  if (setequal(colnames(integrated), cells_before)) {
    integrated <- integrated[, cells_before]
    message("[order] restored input cell order after integration")
  } else {
    stop("integration changed the cell identities, not just their order",
         call. = FALSE)
  }
}

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
message(sprintf("[write] %s", opt$output))
saveRDS(integrated, opt$output, compress = FALSE)
message(sprintf("[done] %.2f GB on disk", file.size(opt$output) / 1024^3))
