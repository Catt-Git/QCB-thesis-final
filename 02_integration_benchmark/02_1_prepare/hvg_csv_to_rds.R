#!/usr/bin/env Rscript

# 02_1_prepare: convert the 01_5 HVG CSV to an .rds character vector, read by
# Seurat CCA/RPCA in 02_2 as anchor.features.
#
# Usage:
#   export DATA_DIR=<data root>
#   Rscript hvg_csv_to_rds.R \
#       -i $DATA_DIR/shiao_hvg_2k_unintegrated_list.csv \
#       -o $DATA_DIR/shiao_hvg_2k_unintegrated_list.rds

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "input HVG CSV (one gene symbol per line, no header)"),
  make_option(c("-o", "--output"), type = "character",
              help = "output .rds"),
  make_option("--n-expected", type = "integer", default = 2000L,
              help = "number of HVGs expected, as a guard [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("both -i/--input and -o/--output are required", call. = FALSE)
}
if (!file.exists(opt$input)) {
  stop(sprintf("input not found: %s", opt$input), call. = FALSE)
}

hvg <- read.csv(opt$input, header = FALSE, stringsAsFactors = FALSE)[[1]]
hvg <- as.character(hvg)

n_expected <- opt[["n-expected"]]
stopifnot(
  "HVG list is empty"            = length(hvg) > 0,
  "unexpected HVG count"         = length(hvg) == n_expected,
  "HVG list has duplicates"      = !any(duplicated(hvg)),
  "HVG list has missing symbols" = !any(is.na(hvg) | hvg == "")
)

saveRDS(hvg, opt$output)
message(sprintf("[done] wrote %d HVG symbols to %s", length(hvg), opt$output))
