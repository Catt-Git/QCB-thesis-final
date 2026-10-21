"""
02_1_prepare: per-batch z-scoring of the HVG object (scaled twin of shiao_hvg_2k.h5ad).

Each cohort is centred and scaled independently, as scIB does. Reimplements
scib.preprocessing.scale_batch, which breaks on anndata >= 0.13 (removed concatenate).

PCA, graph and UMAP from 01_5 describe the unscaled matrix and are dropped.
PCA is recomputed (50 comps, arpack) because the metrics reuse obsm['X_pca'];
the UMAP is recomputed with 01_5's parameters for the 02_3 QC panels; the graph is not kept.

Do not re-run on an existing scaled file that finished integrations already used:
it rewrites .X under them.

Input : $DATA_DIR/shiao_hvg_2k.h5ad         (619,693 x 2,000, sparse)
Output: $DATA_DIR/shiao_hvg_2k_scaled.h5ad  (dense by construction, 1.1 GB)

Usage:
  export DATA_DIR=<data root>
  python3 scale_batch.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import scanpy as sc

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
)
from h5ad_compat import write_h5ad_compat  # noqa: E402

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "shiao_hvg_2k.h5ad")
OUT_PATH = os.path.join(DATA_DIR, "shiao_hvg_2k_scaled.h5ad")

BATCH_KEY = "cohort"
N_PCS = 50
SEED = 0

# Same as uns['neighbors']['params'] of the unscaled object.
N_NEIGHBORS = 15
UMAP_METRIC = "euclidean"

# scanpy scales to a sample std (ddof=1); numpy's default ddof=0 would fail the check.
STD_DDOF = 1
MEAN_TOL = 1e-3
STD_TOL = 1e-2


print(f"Loading {IN_PATH} ...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

assert BATCH_KEY in adata.obs, f"missing batch key {BATCH_KEY!r}"
assert "counts" in adata.layers, "expected raw counts in .layers['counts']"

obs_names_before = adata.obs_names.to_numpy().copy()
var_names_before = adata.var_names.to_numpy().copy()

# Per-batch scaling

batches = adata.obs[BATCH_KEY].cat.categories
print(f"\nScaling {len(batches)} batches of {BATCH_KEY!r} ...", flush=True)

# Results placed by row position, so cell order is preserved.
scaled = np.empty((adata.n_obs, adata.n_vars), dtype=np.float32)
written = np.zeros(adata.n_obs, dtype=bool)

batch_values = adata.obs[BATCH_KEY].to_numpy()
for i, batch in enumerate(batches, start=1):
    rows = np.flatnonzero(batch_values == batch)
    if rows.size == 0:
        raise ValueError(f"batch {batch!r} is an empty category; drop it upstream")

    block = adata.X[rows]
    block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
    block = block.astype(np.float32, copy=False)

    # Same defaults as scIB's scale_batch; zero-variance genes become zeros.
    scaled[rows] = sc.pp.scale(block, zero_center=True, max_value=None, copy=True)
    written[rows] = True

    print(
        f"  [{i:2d}/{len(batches)}] {str(batch):<20} {rows.size:>7,} cells",
        flush=True,
    )

assert written.all(), (
    f"{int((~written).sum())} cells were never assigned to a batch; "
    f"check for NaN in obs[{BATCH_KEY!r}]"
)

adata.X = scaled
del scaled

# Post scaling checks

assert np.array_equal(adata.obs_names.to_numpy(), obs_names_before), (
    "cell order changed during scaling"
)
assert np.array_equal(adata.var_names.to_numpy(), var_names_before), (
    "gene order changed during scaling"
)
assert np.isfinite(adata.X).all(), "scaling produced non-finite values"

# Check per batch: a global zero mean would also pass for whole-matrix scaling.
print("\nVerifying per-batch moments ...", flush=True)
for batch in batches:
    rows = np.flatnonzero(batch_values == batch)
    block = adata.X[rows]
    mean = np.abs(block.mean(axis=0)).max()
    # Genes constant within a batch are all-zero, so excluded.
    std = block.std(axis=0, ddof=STD_DDOF)
    varying = std > 0
    std_error = np.abs(std[varying] - 1).max() if varying.any() else 0.0
    assert mean < MEAN_TOL, f"batch {batch!r}: max |mean| = {mean:.2e}"
    assert std_error < STD_TOL, f"batch {batch!r}: max |std - 1| = {std_error:.2e}"
print(f"  all {len(batches)} batches: |mean| < {MEAN_TOL}, |std - 1| < {STD_TOL}", flush=True)

# Drop everything computed on the unscaled matrix before rebuilding, so a failed
# rebuild leaves nothing stale behind.
for key in ["X_umap"]:
    if key in adata.obsm:
        del adata.obsm[key]
        print(f"\nDropped obsm[{key!r}] (computed on the unscaled matrix)", flush=True)
for key in ["neighbors", "umap"]:
    if key in adata.uns:
        del adata.uns[key]
        print(f"Dropped uns[{key!r}]", flush=True)
for key in list(adata.obsp):
    del adata.obsp[key]
    print(f"Dropped obsp[{key!r}]", flush=True)

print(f"\nRecomputing PCA ({N_PCS} components, arpack) on the scaled matrix ...", flush=True)
# mask_var=None: use all 2,000 genes explicitly, not via var['highly_variable'].
sc.pp.pca(adata, n_comps=N_PCS, svd_solver="arpack", random_state=SEED, mask_var=None)

assert "pca" in adata.uns and "variance" in adata.uns["pca"], "PCA metadata not written"
assert adata.obsm["X_pca"].shape == (adata.n_obs, N_PCS), "unexpected X_pca shape"

# The unintegrated UMAP of this object, for 02_3's QC panels.

print(f"\nComputing the unintegrated UMAP (n_neighbors={N_NEIGHBORS}, "
      f"metric={UMAP_METRIC}, on X_pca) ...", flush=True)
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, metric=UMAP_METRIC,
                use_rep="X_pca", random_state=SEED)
sc.tl.umap(adata, random_state=SEED)

assert adata.obsm["X_umap"].shape == (adata.n_obs, 2), "unexpected X_umap shape"
assert np.isfinite(adata.obsm["X_umap"]).all(), "UMAP produced non-finite coordinates"

# Keep the layout, drop the graph: unused downstream and hundreds of MB.
adata.uns["unintegrated_umap"] = {
    "computed_on": "X_pca",
    "n_neighbors": N_NEIGHBORS,
    "metric": UMAP_METRIC,
    "random_state": SEED,
    "note": ("unintegrated layout of the per-batch scaled matrix; the neighbour "
             "graph behind it is not persisted"),
}
for key in ["neighbors", "umap"]:
    adata.uns.pop(key, None)
for key in list(adata.obsp):
    del adata.obsp[key]
print("  layout kept in obsm['X_umap']; neighbour graph discarded", flush=True)

# Save the scaled object. The output is dense by construction (much larger).

print(f"\nWriting {OUT_PATH} ...", flush=True)
print(f"  .X is dense: {adata.n_obs:,} x {adata.n_vars:,} float32 "
      f"= {adata.X.nbytes / 1024 ** 3:.2f} GB in memory", flush=True)
# Scaled floats barely compress, higher levels only cost time.
write_h5ad_compat(adata, OUT_PATH, compression="gzip", compression_opts=1)
print(
    f"Done: {adata.n_obs:,} x {adata.n_vars:,}, "
    f"{os.path.getsize(OUT_PATH) / 1024 ** 3:.2f} GB on disk",
    flush=True,
)
