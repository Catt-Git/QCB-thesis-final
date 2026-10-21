"""
01_4 celltypist_annotation: cell type annotation with CellTypist

CellTypist requires CP10K+log1p input, so it runs on a temporary copy normalized
from raw counts. Only the predicted labels are kept; the scran .X is untouched.

Input : $DATA_DIR/all_samples_combined_scrublet_norm_cc.h5ad
        Output of 01_4_cell_cycle. .X = scran log-normalized; .layers['counts'] = raw
        counts; gene symbols as var_names.
Output: $DATA_DIR/all_samples_combined_scrublet_norm_cc_annotated.h5ad
        Adds .obs['cell_type'] (majority-voting label, used downstream) and
        .obs['celltypist_predicted'] (raw per-cell prediction, kept for QC).
        .X and .layers['counts'] unchanged.

Model: Cells_Adult_Breast.pkl (Kumar et al. 2023 adult breast atlas), loaded from
       $DATA_DIR, majority_voting=True.

Usage:
  export DATA_DIR=<data root>
  python3 celltypist_annotation.py
"""

from __future__ import annotations
import os
import numpy as np
import scanpy as sc
import anndata as ad
import celltypist
from celltypist import models

sc.settings.verbosity = 1

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm_cc.h5ad")
OUT_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm_cc_annotated.h5ad")

# Model loaded from the .pkl in DATA_DIR, not from the CellTypist cache.
MODEL_PATH = os.path.join(DATA_DIR, "Cells_Adult_Breast.pkl")

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

assert "counts" in adata.layers, "Expected raw counts in .layers['counts']"

print(f"Loading CellTypist model from: {MODEL_PATH}", flush=True)
assert os.path.exists(MODEL_PATH), f"Model file not found: {MODEL_PATH}"
model = models.Model.load(model=MODEL_PATH)

# Temporary CP10K+log1p copy from raw counts (the normalization CellTypist was trained on).
print("Building temporary CP10K+log1p matrix for CellTypist input...", flush=True)
adata_ct = adata.copy()
adata_ct.X = adata_ct.layers["counts"].copy()
sc.pp.normalize_total(adata_ct, target_sum=1e4)
sc.pp.log1p(adata_ct)

print("Annotating with CellTypist (majority_voting=True, CPU)...", flush=True)
# majority_voting refines per-cell predictions over over-clustered neighborhoods.
predictions = celltypist.annotate(
    adata_ct,
    model=model,
    majority_voting=True,
    use_GPU=False,
)

# cell_type = majority-voting label (used downstream); celltypist_predicted = per-cell call, for QC.
adata.obs["cell_type"] = predictions.predicted_labels["majority_voting"].values
adata.obs["celltypist_predicted"] = predictions.predicted_labels["predicted_labels"].values

del adata_ct  # drop the temporary normalized copy

print("cell_type label distribution:", flush=True)
print(adata.obs["cell_type"].value_counts(), flush=True)

# Post-condition checks: labels added, main matrices untouched.
assert "cell_type" in adata.obs
assert "counts" in adata.layers, "raw counts layer lost"
assert "celltypist" not in adata.layers, "temporary celltypist layer must not be persisted"

print("Saving...", flush=True)
adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)