"""
03_2 scran normalization of the malignant subset (mirrors 01_3/scran_norm.py).

Removing the non-malignant cells changes the scran pools, so the subset is re-normalized
from raw counts with the parameters of 01_3.

Input : $DATA_DIR/03_tum/shiao_tum_raw.h5ad
        Output of subset_and_qc.ipynb. .X = .layers['counts'] = raw counts; zero-count
        genes already removed.
Output: $DATA_DIR/03_tum/shiao_tum_norm.h5ad
        .X = scran log-normalized, float32; .layers['counts'] unchanged;
        .obs['size_factors'] added; .raw cleared.

Usage:
  export DATA_DIR=<data root>
  ./subsetting_all.sh norm       # step 1 of 3
  python3 scran_norm_tum.py      # standalone, no logging or resume
"""

from __future__ import annotations
import os
import sys
import numpy as np
import scanpy as sc
import scipy.sparse as _sp

import anndata2ri
from rpy2.robjects import conversion, default_converter

# Recent anndata2ri dropped activate/deactivate, which scib still calls.
if not hasattr(anndata2ri, "activate"):
    _a2r_converter = default_converter + anndata2ri.converter

    def _activate():
        conversion.set_conversion(_a2r_converter)

    def _deactivate():
        conversion.set_conversion(default_converter)

    anndata2ri.activate = _activate
    anndata2ri.deactivate = _deactivate

# scib's normalize leaves a COO .X (sparse .multiply), which anndata >=0.11 rejects:
# coerce the result back to CSR, as in 01_3.
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
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
import seaborn as sns

# cell_set lives in the phase's utils/.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import cell_set as C

# scib's normalize has no random_state, so the global RNGs of its leiden preclustering are seeded.
random.seed(C.SEED)
np.random.seed(C.SEED)
sc.settings.verbosity = 1

C.banner("03_2 scran normalization")
IN_PATH = C.path("_raw.h5ad")
OUT_PATH = C.path("_norm.h5ad")

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

assert "counts" in adata.layers, "Expected raw counts in .layers['counts']"

# Guard against re-normalizing log-normalized values: .X must be raw counts here.
x_sample = adata.X[:100].data if hasattr(adata.X[:100], "data") else np.asarray(adata.X[:100]).ravel()
assert np.allclose(x_sample, np.round(x_sample)), (
    ".X does not look like raw integer counts; subset_and_qc.ipynb must restore "
    ".X = .layers['counts'] before this step"
)

n_zero_genes = int(np.asarray((adata.X.sum(axis=0) == 0)).sum())
assert n_zero_genes == 0, f"{n_zero_genes} zero-count genes present; run sc.pp.filter_genes(min_cells=1) first"

assert "size_factors" not in adata.obs, (
    "obs['size_factors'] is still present; it was estimated on a larger population "
    "and must be dropped by subset_and_qc.ipynb"
)

print("Normalizing data (scran via scib)...", flush=True)
# quickCluster -> computeSumFactors -> logNormCounts, as in 01_3.
scib.preprocessing.normalize(
    adata,
    min_mean=0.1,
    log=True,
    precluster=True,
    cluster_method="leiden",
    sparsify=False,
)

# scib leaves an unused .raw copy and a float64 .X.
adata.raw = None
if adata.X.dtype != np.float32:
    print(f"Downcasting .X from {adata.X.dtype} to float32", flush=True)
    adata.X = adata.X.astype(np.float32)

assert "size_factors" in adata.obs, "scib did not write size_factors"
assert np.all(np.isfinite(adata.obs["size_factors"].values)), "non-finite size factors present"
assert "counts" in adata.layers, "raw counts layer lost during normalization"
assert adata.raw is None, ".raw should not be set"
assert adata.X.dtype == np.float32, ".X should be float32"

print("Saving normalized data...", flush=True)
adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)

print("Generating normalization diagnostics...", flush=True)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PHASE_DIR = os.path.dirname(SCRIPT_DIR)  # 03_drvi_tumoral_epi/
FIG_DIR = os.environ.get("FIG_DIR", os.path.join(PHASE_DIR, "figures", "03_2_normalization"))
os.makedirs(FIG_DIR, exist_ok=True)

sf = adata.obs["size_factors"].to_numpy()
n_nonpos = int((sf <= 0).sum())
print(f"size_factors: min={sf.min():.4f}, max={sf.max():.4f}, n(<=0)={n_nonpos}", flush=True)
if n_nonpos > 0:
    print(f"WARNING: {n_nonpos} cells with non-positive size factors", flush=True)

rng = np.random.default_rng(C.SEED)
n_plot = min(30000, adata.n_obs)
idx = rng.choice(adata.n_obs, size=n_plot, replace=False)
sub = adata.obs.iloc[idx]

# Coloured by treatment: cell_type has a single level on the malignant subset.
suffix = C.compartment()
fig, ax = plt.subplots(figsize=(6, 6))
hue = C.TREATMENT_KEY if C.TREATMENT_KEY in sub.columns else None
sns.scatterplot(data=sub, x="total_counts", y="size_factors",
                hue=hue, s=4, alpha=0.3, linewidth=0, ax=ax)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Total raw counts (library size)")
ax.set_ylabel("scran size factors")
ax.set_title(f"Size factors vs library size ({suffix})")
ax.grid(True, which="both", ls="-", alpha=0.2)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, f"scran_size_factors_vs_library_size_{suffix}.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4))
sns.histplot(sf, bins=50, kde=False, ax=ax)
ax.set_xlabel("scran size factors"); ax.set_ylabel("Number of cells")
ax.set_title(f"Distribution of scran size factors ({suffix})")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, f"scran_size_factors_distribution_{suffix}.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved diagnostics to {FIG_DIR}", flush=True)
