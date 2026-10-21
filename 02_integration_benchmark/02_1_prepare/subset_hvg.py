"""
02_1_prepare: restrict the unintegrated object to the 2,000 HVGs, once for all methods.

- layers['counts'] must stay raw integers: scVI, scANVI and DRVI train on it.
- X_pca and uns['pca'] from 01_5 (computed on these same genes) are kept, so
  scib's pcr() reuses them as the shared baseline.
Genes keep the source object's order (boolean mask, not CSV order).

Input : $DATA_DIR/shiao.h5ad                     (from 01_5, 619,693 x 30,869)
        $DATA_DIR/shiao_hvg_2k_unintegrated_list.csv  (from 01_5, 2,000 symbols)
Output: $DATA_DIR/shiao_hvg_2k.h5ad              (619,693 x 2,000, 778 MB, sparse)

Usage:
  export DATA_DIR=<data root>
  python3 subset_hvg.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
)
from h5ad_compat import write_h5ad_compat  # noqa: E402

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "shiao.h5ad")
HVG_PATH = os.path.join(DATA_DIR, "shiao_hvg_2k_unintegrated_list.csv")
OUT_PATH = os.path.join(DATA_DIR, "shiao_hvg_2k.h5ad")

N_HVGS = 2000
BATCH_KEY = "cohort"
LABEL_KEY = "cell_type"


def sparse_values(matrix) -> np.ndarray:
    """The stored values of a matrix, sparse or dense, as a flat array."""
    return matrix.data if hasattr(matrix, "data") else np.asarray(matrix).ravel()


def describe(matrix, name: str) -> None:
    values = sparse_values(matrix)
    nnz = values.size
    total = matrix.shape[0] * matrix.shape[1]
    print(
        f"  {name:<8} dtype={matrix.dtype}  nnz={nnz:,}  "
        f"density={nnz / total:.3%}  min={values.min():.4g}  max={values.max():.4g}",
        flush=True,
    )


print(f"Loading {IN_PATH} ...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

# Input check

assert "counts" in adata.layers, "expected raw counts in .layers['counts']"
assert BATCH_KEY in adata.obs, f"missing batch key {BATCH_KEY!r}"
assert LABEL_KEY in adata.obs, f"missing label key {LABEL_KEY!r}"
assert "highly_variable" in adata.var, "missing var['highly_variable'] from 01_5"

hvg_csv = pd.read_csv(HVG_PATH, header=None)[0].astype(str).tolist()
assert len(hvg_csv) == N_HVGS, f"{HVG_PATH} has {len(hvg_csv)} rows, expected {N_HVGS}"

# Compared as sets, so the check does not depend on the order of the CSV.
flagged = set(adata.var_names[adata.var["highly_variable"].to_numpy()])
missing = set(hvg_csv) - flagged
extra = flagged - set(hvg_csv)
assert not missing and not extra, (
    "the HVG csv and var['highly_variable'] disagree: "
    f"{len(missing)} only in the csv, {len(extra)} only in the flag. "
    "They come from the same 01_5 run and must be identical."
)

# Subset

mask = np.asarray(adata.var_names.isin(hvg_csv))
assert mask.sum() == N_HVGS, f"mask selects {mask.sum()} genes, expected {N_HVGS}"

n_obs_before = adata.n_obs
obs_names_before = adata.obs_names.to_numpy().copy()

print(f"Subsetting to {N_HVGS} HVGs ...", flush=True)
adata = adata[:, mask].copy()

# Post-subset checks

assert adata.n_vars == N_HVGS, f"got {adata.n_vars} genes after subsetting"
assert adata.n_obs == n_obs_before, "cells were lost during the gene subset"
assert np.array_equal(adata.obs_names.to_numpy(), obs_names_before), (
    "cell order changed during the gene subset"
)
assert bool(adata.var["highly_variable"].all()), "a non-HVG survived the subset"

# Raw counts must still be raw counts.
counts = adata.layers["counts"]
counts_values = sparse_values(counts)
assert np.allclose(counts_values, np.floor(counts_values)), (
    "layers['counts'] contains non-integer values: it is no longer raw counts. "
    "scVI, scANVI and DRVI read this layer directly, so stop here rather than "
    "train on normalised data."
)
assert counts_values.min() >= 0, "layers['counts'] contains negative values"

# And .X must not be raw counts (catches a swap with the counts layer upstream).
x_values = sparse_values(adata.X)
assert not np.allclose(x_values, np.floor(x_values)), (
    ".X looks like integer counts, but it should be scran log-normalised "
    "expression. Check whether .X and layers['counts'] were swapped upstream."
)

# PCA baseline reused by scib's pcr().
assert "X_pca" in adata.obsm, "obsm['X_pca'] missing, the metrics baseline needs it"
assert "pca" in adata.uns and "variance" in adata.uns["pca"], (
    "uns['pca']['variance'] missing: scib would silently recompute the PCA in "
    "every metrics job instead of reusing this one"
)
assert adata.obsm["X_pca"].shape[0] == adata.n_obs, "X_pca has the wrong number of rows"

print("\nContents of the subset object:", flush=True)
describe(adata.X, ".X")
describe(adata.layers["counts"], "counts")
print(f"  obsm     {list(adata.obsm)}", flush=True)
print(f"  varm     {list(adata.varm)}", flush=True)
print(f"  obsp     {list(adata.obsp)}", flush=True)
print(f"  X_pca    {adata.obsm['X_pca'].shape}", flush=True)
print(f"  batches  {adata.obs[BATCH_KEY].nunique()} in {BATCH_KEY!r}", flush=True)
print(f"  labels   {adata.obs[LABEL_KEY].nunique()} in {LABEL_KEY!r}", flush=True)

# Save the subset object

print(f"\nWriting {OUT_PATH} ...", flush=True)
write_h5ad_compat(adata, OUT_PATH, compression="gzip")
print(
    f"Done: {adata.n_obs:,} x {adata.n_vars:,}, "
    f"{os.path.getsize(OUT_PATH) / 1024 ** 3:.2f} GB on disk",
    flush=True,
)
