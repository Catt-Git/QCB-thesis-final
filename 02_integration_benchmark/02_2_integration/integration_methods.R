# 02_2 integration: R method implementations, adapted from Luecken et al. 2022
# (Nat Methods, https://doi.org/10.1038/s41592-021-01336-8).
#
#   fastmnn      (batchelor)  -> corrected RNA assay (full) + fastmnn reduction (embed)
#   seurat_cca   (Seurat)     -> corrected 'integrated' assay (full)
#   seurat_rpca  (Seurat)     -> corrected 'integrated' assay (full)
#
# Differences from the original:
#   - Seurat 5.5.1: no do.cpp in IntegrateData; RPCA gets the HVG vector, not 2000.
#   - Legacy FindIntegrationAnchors + IntegrateData kept, so Seurat stays a full method.
#   - CCA/RPCA are reference-based (Patient53/16/43 by default): pairwise anchoring
#     of 34 batches is not feasible at 620k cells.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

# Reference patient names -> positions in the SplitObject list.
.reference_indices <- function(batch_list, reference_batches) {
  idx <- match(reference_batches, names(batch_list))
  if (any(is.na(idx))) {
    missing <- reference_batches[is.na(idx)]
    stop(sprintf("reference batch(es) not found in the data: %s\navailable: %s",
                 paste(missing, collapse = ", "),
                 paste(names(batch_list), collapse = ", ")), call. = FALSE)
  }
  idx
}

# Seurat CCA; corrected matrix in assay 'integrated'.

runSeuratCCA <- function(sobj, batch, hvg, reference_batches,
                         dims = 1:30, k.weight = 100) {
  batch_list <- SplitObject(sobj, split.by = batch)
  ref_idx <- .reference_indices(batch_list, reference_batches)
  message(sprintf("[cca] %d batches, reference: %s",
                  length(batch_list), paste(reference_batches, collapse = ", ")))

  anchors <- FindIntegrationAnchors(
    object.list     = batch_list,
    anchor.features = hvg,
    reference       = ref_idx,
    reduction       = "cca",
    scale           = TRUE,
    l2.norm         = TRUE,
    dims            = dims,
    k.anchor        = 5,
    k.filter        = 200,
    k.score         = 30,
    max.features    = 200,
    eps             = 0
  )
  IntegrateData(
    anchorset        = anchors,
    new.assay.name   = "integrated",
    dims             = dims,
    k.weight         = k.weight,
    sd.weight        = 1,
    preserve.order   = FALSE,
    eps              = 0,
    verbose          = TRUE
  )
}

# Seurat RPCA; each batch is scaled and PCA'd on the HVGs before anchoring.

runSeuratRPCA <- function(sobj, batch, hvg, reference_batches,
                          dims = 1:30, k.weight = 100) {
  batch_list <- SplitObject(sobj, split.by = batch)
  ref_idx <- .reference_indices(batch_list, reference_batches)
  message(sprintf("[rpca] %d batches, reference: %s",
                  length(batch_list), paste(reference_batches, collapse = ", ")))

  # features = hvg must be the vector of HVG symbols, not a count.
  batch_list <- lapply(batch_list, function(x) {
    x <- ScaleData(x, features = hvg, verbose = FALSE)
    x <- RunPCA(x, features = hvg, verbose = FALSE)
    x
  })

  anchors <- FindIntegrationAnchors(
    object.list     = batch_list,
    anchor.features = hvg,
    reference       = ref_idx,
    reduction       = "rpca",
    scale           = TRUE,
    l2.norm         = TRUE,
    dims            = dims,
    k.anchor        = 5,
    k.filter        = 200,
    k.score         = 30,
    max.features    = 200,
    eps             = 0
  )
  IntegrateData(
    anchorset        = anchors,
    new.assay.name   = "integrated",
    dims             = dims,
    k.weight         = k.weight,
    sd.weight        = 1,
    preserve.order   = FALSE,
    eps              = 0,
    verbose          = TRUE
  )
}

# fastMNN on the 'data' layer: corrected expression in RNA (full), embedding in
# the 'fastmnn' reduction (embed).

runFastMNN <- function(sobj, batch) {
  suppressPackageStartupMessages({
    library(batchelor)
    library(Matrix)
    library(SingleCellExperiment)
  })

  expr <- GetAssayData(sobj, assay = DefaultAssay(sobj), layer = "data")
  sce <- fastMNN(expr, batch = sobj@meta.data[[batch]])

  # Data slot only: a counts copy of the dense matrix would double the .rds.
  corrected <- as(assay(sce, "reconstructed"), "CsparseMatrix")
  sobj[["RNA"]] <- CreateAssayObject(data = corrected)
  DefaultAssay(sobj) <- "RNA"

  # Corrected embedding -> reduction (embed output).
  emb <- reducedDim(sce, "corrected")
  rownames(emb) <- colnames(sobj)
  sobj[["fastmnn"]] <- CreateDimReducObject(
    embeddings = emb, key = "fastmnn_", assay = "RNA"
  )
  sobj
}

# Registry consumed by run_integration.R.
INTEGRATION_METHODS <- list(
  fastmnn     = function(sobj, batch, hvg, reference_batches)
    runFastMNN(sobj, batch),
  seurat_cca  = function(sobj, batch, hvg, reference_batches)
    runSeuratCCA(sobj, batch, hvg, reference_batches),
  seurat_rpca = function(sobj, batch, hvg, reference_batches)
    runSeuratRPCA(sobj, batch, hvg, reference_batches)
)
