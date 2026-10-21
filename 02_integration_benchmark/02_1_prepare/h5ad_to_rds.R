#!/usr/bin/env Rscript

# 02_1_prepare: convert an .h5ad benchmark input to a Seurat v3 .rds for the R
# methods (fastMNN, Seurat CCA/RPCA) in 02_2.
#
# Read with zellkonverter's native R reader (no Python). The object is forced to
# v3 because the legacy FindIntegrationAnchors + IntegrateData path needs it.
#   counts slot <- layers['counts']  (raw, sparse; Seurat v3 needs non-negative counts)
#   data slot   <- .X                (log-normalised, or z-scored for the scaled input)
# Sparse/dense storage is kept as read: the scaled .X is genuinely dense.
#
# Usage:
#   export DATA_DIR=<data root>
#   Rscript h5ad_to_rds.R -i $DATA_DIR/shiao_hvg_2k.h5ad        -o $DATA_DIR/shiao_hvg_2k.rds
#   Rscript h5ad_to_rds.R -i $DATA_DIR/shiao_hvg_2k_scaled.h5ad -o $DATA_DIR/shiao_hvg_2k_scaled.rds

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "input .h5ad"),
  make_option(c("-o", "--output"), type = "character",
              help = "output .rds"),
  make_option("--x-name", type = "character", default = "logcounts",
              help = "name to give the .X assay after reading [default %default]"),
  make_option("--counts-layer", type = "character", default = "counts",
              help = "layer holding raw counts for the Seurat counts slot [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("both -i/--input and -o/--output are required", call. = FALSE)
}
if (!file.exists(opt$input)) {
  stop(sprintf("input not found: %s", opt$input), call. = FALSE)
}

x_name <- opt[["x-name"]]
counts_layer <- opt[["counts-layer"]]

# Read

message(sprintf("[read] %s", opt$input))
t0 <- Sys.time()
sce <- readH5AD(opt$input, reader = "R")
message(sprintf("[read] %.1f s, %d genes x %d cells, assays: %s",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                nrow(sce), ncol(sce),
                paste(assayNames(sce), collapse = ", ")))

# The R reader names the .X assay "X" regardless of X_name; rename it.
if (!"X" %in% assayNames(sce)) {
  stop(sprintf("expected an assay named 'X' from the R reader, found: %s",
               paste(assayNames(sce), collapse = ", ")), call. = FALSE)
}
names(assays(sce))[names(assays(sce)) == "X"] <- x_name

has_counts <- counts_layer %in% assayNames(sce)
if (!has_counts) {
  message(sprintf("[warn] no '%s' assay; using %s as the counts slot too",
                  counts_layer, x_name))
}

# Build V3 Seurat object

options(Seurat.object.assay.version = "v3")

# Sparse -> dgCMatrix, dense stays dense (scaled .X would be larger as sparse).
as_seurat_matrix <- function(m) {
  if (is(m, "sparseMatrix")) as(m, "CsparseMatrix") else as.matrix(m)
}

data_mat <- as_seurat_matrix(assay(sce, x_name))
counts_mat <- if (has_counts) as_seurat_matrix(assay(sce, counts_layer)) else data_mat
message(sprintf("[build] data slot is %s, counts slot is %s",
                if (is(data_mat, "sparseMatrix")) "sparse" else "dense",
                if (is(counts_mat, "sparseMatrix")) "sparse" else "dense"))

message("[build] assembling v3 Seurat object")
so <- CreateSeuratObject(
  counts    = counts_mat,
  meta.data = as.data.frame(colData(sce)),
  assay     = "RNA"
)
so <- SetAssayData(so, layer = "data", new.data = data_mat)

# Carry over X_pca / X_umap so downstream code need not recompute them.
for (dr in reducedDimNames(sce)) {
  emb <- reducedDim(sce, dr)
  rownames(emb) <- colnames(so)
  key <- paste0(gsub("[^A-Za-z0-9]", "", tolower(dr)), "_")
  so[[dr]] <- CreateDimReducObject(embeddings = emb, key = key, assay = "RNA")
}

# Check the Seurat object is consistent with the SCE it came from
stopifnot(
  "assay is not v3"          = inherits(so[["RNA"]], "Assay"),
  "cell order changed"       = identical(colnames(so), colnames(sce)),
  "gene order changed"       = identical(rownames(so), rownames(sce)),
  "data slot lost its shape" = identical(dim(GetAssayData(so, layer = "data")),
                                         dim(assay(sce, x_name)))
)
if (has_counts) {
  cnts <- GetAssayData(so, layer = "counts")
  stopifnot("counts slot holds non-integers" = all(cnts@x == floor(cnts@x)),
            "counts slot holds negatives"     = min(cnts@x) >= 0)
}

# Write the Seurat object to .rds
# compress=FALSE: gzip on the ~10 GB dense scaled object is slow for little gain.
message(sprintf("[write] %s", opt$output))
t1 <- Sys.time()
saveRDS(so, opt$output, compress = FALSE)
message(sprintf("[write] %.1f s, %.2f GB on disk",
                as.numeric(difftime(Sys.time(), t1, units = "secs")),
                file.size(opt$output) / 1024^3))
message("[done]")
