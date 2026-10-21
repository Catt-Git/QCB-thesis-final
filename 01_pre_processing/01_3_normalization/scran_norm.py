"""
01_3_normalization: scran normalization via scib wrapper

Usage:
  export DATA_DIR=<data root>
  python3 scran_norm.py

Input : $DATA_DIR/all_samples_combined_scrublet.h5ad
        AnnData after 01_2 (QC + Scrublet + filters).
        .X = raw counts (int); .layers['counts'] = identical raw counts.
        Genes detected in fewer than 3 cells already removed upstream (min_cells=3).
Output: $DATA_DIR/all_samples_combined_scrublet_norm.h5ad
        .X = scran log-normalized, float32; .layers['counts'] = raw counts (unchanged);
        .obs['size_factors'] added by scib. .raw is explicitly cleared (scib sets
        adata.raw = adata, a full unused duplicate of .X -- see comment below).
"""

from __future__ import annotations
import os
import numpy as np
import scanpy as sc
import anndata as ad
import scipy.sparse as _sp

import anndata2ri
from rpy2.robjects import conversion, default_converter

if not hasattr(anndata2ri, "activate"):
    _a2r_converter = default_converter + anndata2ri.converter

    def _activate():
        conversion.set_conversion(_a2r_converter)

    def _deactivate():
        conversion.set_conversion(default_converter)

    anndata2ri.activate = _activate
    anndata2ri.deactivate = _deactivate

# Patch for scib 1.1.7: its normalize step uses sparse .multiply(), which returns a COO
# matrix, and anndata >= 0.11 only accepts CSR/CSC in .X, so scib crashes.
# Convert multiply's output back to CSR.

for _cls in (_sp.csr_matrix, _sp.csc_matrix):
    if not getattr(_cls.multiply, "_csr_coerced", False):
        _orig_multiply = _cls.multiply

        def _multiply_tocsr(self, other, _orig=_orig_multiply):
            result = _orig(self, other)
            return result.tocsr() if _sp.issparse(result) else result

        _multiply_tocsr._csr_coerced = True
        _cls.multiply = _multiply_tocsr

import scib
import random
import matplotlib
matplotlib.use("Agg")  # headless backend for SLURM/no-display runs
import matplotlib.pyplot as plt
import seaborn as sns

SEED = 0
random.seed(SEED)
np.random.seed(SEED)
sc.settings.verbosity = 1

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet.h5ad")
OUT_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm.h5ad")

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

# Sanity checks on the input object (do not silently proceed if violated)
assert "counts" in adata.layers, "Expected raw counts in .layers['counts']"
# .X must be raw integer counts at this point (scran normalizes from raw)
x_sample = adata.X[:100].data if hasattr(adata.X[:100], "data") else np.asarray(adata.X[:100]).ravel()
assert np.allclose(x_sample, np.round(x_sample)), ".X does not look like raw integer counts"

# Guard against re-running the zero-count gene filter: it was applied in 01_2.
# scran's computeSumFactors fails on all-zero genes, so we assert rather than refilter.
n_zero_genes = int(np.asarray((adata.X.sum(axis=0) == 0)).sum())
assert n_zero_genes == 0, f"{n_zero_genes} zero-count genes present; run the 01_2 gene filter (min_cells=3) first"

print("Normalizing data (scran via scib)...", flush=True)
# scib scran wrapper: quickCluster -> computeSumFactors -> logNormCounts.
# Raw counts remain in .layers['counts']; size factors written to .obs['size_factors'].
scib.preprocessing.normalize(
    adata,
    min_mean=0.1, # computeSumFactors min.mean; standard scib default for droplet data
    log=True, # return log1p-transformed normalized expression in .X
    precluster=True, # run quickCluster (leiden) before size-factor estimation
    cluster_method="leiden", # leiden preclustering (louvain deprecated)
    sparsify=False, # keep .X storage as-is, no forced conversion
)

# scib leaves .raw as a full copy of .X (unused downstream) and casts .X to float64.
# Drop .raw and cast .X back to float32: the output goes from ~27GB to ~10GB.
adata.raw = None
if adata.X.dtype != np.float32:
    print(f"Downcasting .X from {adata.X.dtype} to float32", flush=True)
    adata.X = adata.X.astype(np.float32)

# Post-condition checks
assert "size_factors" in adata.obs, "scib did not write size_factors"
assert np.all(np.isfinite(adata.obs["size_factors"].values)), "non-finite size factors present"
assert "counts" in adata.layers, "raw counts layer lost during normalization"
assert adata.raw is None, ".raw should not be set (see comment above)"
assert adata.X.dtype == np.float32, ".X should be float32"

print("Saving normalized data...", flush=True)
adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)


print("Generating normalization diagnostics...", flush=True)
# Diagnostic figures go to 01_pre_processing/figures/, anchored to this script's location.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PHASE_DIR = os.path.dirname(SCRIPT_DIR)  # 01_pre_processing/
FIG_DIR = os.environ.get("FIG_DIR", os.path.join(PHASE_DIR, "figures", "01_3_normalization"))
os.makedirs(FIG_DIR, exist_ok=True)

# Guard: scran size factors must be strictly positive (<= 0 breaks the log and
# flags pathological cells that survived QC). Report before plotting.
sf = adata.obs["size_factors"].to_numpy()
n_nonpos = int((sf <= 0).sum())
print(f"size_factors: min={sf.min():.4f}, max={sf.max():.4f}, n(<=0)={n_nonpos}", flush=True)
if n_nonpos > 0:
    print(f"WARNING: {n_nonpos} cells with non-positive size factors", flush=True)

# Subsample for scatter readability (plot only; data untouched)
rng = np.random.default_rng(0)
n_plot = min(30000, adata.n_obs)
idx = rng.choice(adata.n_obs, size=n_plot, replace=False)
sub = adata.obs.iloc[idx]

# First plot:size_factors vs total_counts (log-log), colored by technical fraction.
#    Expectation: monotone positive relation; deviations flag odd cells/clusters.
fig, ax = plt.subplots(figsize=(6, 6))
hue = "dataset_origin" if "dataset_origin" in sub.columns else None
sns.scatterplot(data=sub, x="total_counts", y="size_factors",
                hue=hue, s=4, alpha=0.3, linewidth=0, ax=ax)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Total raw counts (library size)")
ax.set_ylabel("scran size factors")
ax.set_title("Size factors vs library size")
ax.grid(True, which="both", ls="-", alpha=0.2)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "scran_size_factors_vs_library_size.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)

# Second plot: Distribution of size factors: expect unimodal, positive, no extreme tails.
fig, ax = plt.subplots(figsize=(6, 4))
sns.histplot(sf, bins=50, kde=False, ax=ax)
ax.set_xlabel("scran size factors")
ax.set_ylabel("Number of cells")
ax.set_title("Distribution of scran size factors")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "scran_size_factors_distribution.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved diagnostics to {FIG_DIR}", flush=True)