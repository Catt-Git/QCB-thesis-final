"""
03_2 cell cycle scoring of the malignant subset (mirrors 01_4/cell_cycle_score.py).

Re-scored, not inherited: the reference bins of score_genes_cell_cycle depend on the cells
in the object, and this subset is enriched for proliferating cells.

Input : $DATA_DIR/03_tum/shiao_tum_norm.h5ad
        .X = scran log-normalized.
        $CC_GENES
        Tirosh/Regev list (43 S + 54 G2M), one symbol per line.
Output: $DATA_DIR/03_tum/shiao_tum_norm_cc.h5ad
        Adds .obs['S_score'], .obs['G2M_score'], .obs['phase'].

Usage:
  export DATA_DIR=<data root>
  ./subsetting_all.sh cc         # step 2 of 3, sets $CC_GENES

  export CC_GENES=$DATA_DIR/regev_lab_cell_cycle_genes.txt
  python3 cell_cycle_score_tum.py
"""

from __future__ import annotations
import os
import sys
import numpy as np
import scanpy as sc

# cell_set lives in the phase's utils/.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import cell_set as C

sc.settings.verbosity = 1

C.banner("03_2 cell cycle scoring")
IN_PATH = C.path("_norm.h5ad")
OUT_PATH = C.path("_norm_cc.h5ad")
CC_GENES = os.environ["CC_GENES"]

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

# Scoring needs log-normalized .X: raw counts would be all integers.
x_sample = adata.X[:100].data if hasattr(adata.X[:100], "data") else np.asarray(adata.X[:100]).ravel()
x_sample = x_sample[np.isfinite(x_sample)]
assert not np.allclose(x_sample, np.round(x_sample)), \
    ".X looks like raw integer counts; cell cycle scoring expects log-normalized .X"

with open(CC_GENES) as f:
    cell_cycle_genes = [x.strip() for x in f if x.strip()]
assert len(cell_cycle_genes) == 97, \
    f"expected 97 cell cycle genes (43 S + 54 G2M), got {len(cell_cycle_genes)}"
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

# The subset was re-filtered, so some cell cycle genes may be gone: report the survivors.
s_present = [g for g in s_genes if g in adata.var_names]
g2m_present = [g for g in g2m_genes if g in adata.var_names]
print(f"S genes:   {len(s_present)}/{len(s_genes)} present", flush=True)
print(f"G2M genes: {len(g2m_present)}/{len(g2m_genes)} present", flush=True)
missing = sorted(set(s_genes + g2m_genes) - set(s_present) - set(g2m_present))
if missing:
    print(f"Missing cell cycle genes (not in var_names): {missing}", flush=True)

print("Scoring cell cycle...", flush=True)
sc.tl.score_genes_cell_cycle(
    adata,
    s_genes=s_present,
    g2m_genes=g2m_present,
    random_state=C.SEED,
)

print("Cell cycle phase distribution:", flush=True)
print(adata.obs["phase"].value_counts(), flush=True)

# Cycling fraction per pre-CNV label, for reference only.
if C.PRIOR_LABEL_KEY in adata.obs:
    print("\nphase by pre-CNV CellTypist label:", flush=True)
    print(adata.obs.groupby(C.PRIOR_LABEL_KEY, observed=True)["phase"]
          .value_counts(normalize=True).unstack().round(3).to_string(), flush=True)

for col in ("S_score", "G2M_score", "phase"):
    assert col in adata.obs, f"scoring did not add {col}"

print("\nSaving...", flush=True)
adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)
