"""
03_2 feature selection + dimensionality reduction of the malignant subset
(mirrors 01_5/scib_reduce_data.py), and the DRVI input of 03_3.

HVGs are re-selected inside the tumor (the 01_5 ones mostly separate cell types), with
01_5's scib settings; MT- genes are replaced by the next genes of the ranking.

Input : $DATA_DIR/03_tum/shiao_tum_norm_cc.h5ad
        .X = scran log-normalized; .layers['counts'] = raw counts;
        .obs['cohort'] (batch key), 'cell_type', 'cell_type_01_4', 'cnv_status', ...
Output: $DATA_DIR/03_tum/shiao_tum.h5ad
        All genes, .var['highly_variable'] (panel without MT-), PCA, neighbours, UMAP.
        $DATA_DIR/03_tum/shiao_tum_hvg_2k_nomt.h5ad + _hvg_2k_nomt_list.csv
        The same object on the panel only: the DRVI input of 03_3.

Usage:
  export DATA_DIR=<data root>
  ./subsetting_all.sh reduce     # step 3 of 3
  python3 reduce_data_tum.py     # standalone
"""

from __future__ import annotations
import os
import sys
import matplotlib
matplotlib.use("Agg")  # headless backend
import numpy as np
import pandas as pd

# scib's hvg_batch uses np.in1d, removed in numpy>=2.4.
if not hasattr(np, "in1d"):
    np.in1d = np.isin

import anndata2ri
from rpy2.robjects import conversion, default_converter

# Recent anndata2ri dropped activate/deactivate, which scib still calls.
def _activate():
    conversion.set_conversion(default_converter + anndata2ri.converter)

def _deactivate():
    conversion.set_conversion(default_converter)

anndata2ri.activate = _activate
anndata2ri.deactivate = _deactivate

import scanpy as sc
import scib

# cell_set lives in the phase's utils/.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import cell_set as C

# Exact prefix, not a regex: MT1A / MT2A / MTRNR2L12 are nuclear and must not match.
MT_PREFIX = "MT-"

C.banner("03_2 feature selection + reduction")
IN_PATH = C.path("_norm_cc.h5ad")
OUT_PATH = C.path(".h5ad")
HVG_CSV_PATH = C.hvg_path("_list.csv")
HVG_H5AD_PATH = C.hvg_path(".h5ad")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PHASE_DIR = os.path.dirname(SCRIPT_DIR)  # 03_drvi_tumoral_epi/
FIG_DIR = os.environ.get("FIG_DIR", os.path.join(PHASE_DIR, "figures", "03_2_reduce_data"))
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=300, facecolor="white")


def hvg_stats(adata) -> pd.DataFrame:
    """Per-gene statistics of scib's hvg_batch: scanpy's cell_ranger selection per cohort."""
    a = sc.AnnData(X=adata.X, obs=adata.obs[[C.BATCH_KEY]], var=pd.DataFrame(index=adata.var_names))
    sc.pp.highly_variable_genes(a, flavor="cell_ranger", n_top_genes=C.N_HVGS, n_bins=20,
                                batch_key=C.BATCH_KEY)
    return a.var[["dispersions_norm", "highly_variable_nbatches"]]


def ranked_hvgs(stats: pd.DataFrame, n_batches: int, target: int) -> list[str]:
    """hvg_batch's fill order, readable past 2,000 genes."""
    selected: list[str] = []
    for k in range(n_batches, 0, -1):
        tier = stats.loc[stats["highly_variable_nbatches"] == k, "dispersions_norm"]
        selected.extend(tier.sort_values(ascending=False).index[: target - len(selected)])
        if len(selected) >= target:
            break
    if len(selected) < target:
        raise SystemExit(f"only {len(selected)} genes are highly variable in at least one "
                         f"cohort; cannot reach {target}")
    return selected


adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

assert C.BATCH_KEY in adata.obs, f"missing batch key {C.BATCH_KEY!r}"
assert C.LABEL_KEY in adata.obs, f"missing label key {C.LABEL_KEY!r}"
assert "counts" in adata.layers, "expected raw counts in .layers['counts']"

cells_per_batch = adata.obs[C.BATCH_KEY].value_counts()
cells_per_batch = cells_per_batch[cells_per_batch > 0]
print(f"{C.BATCH_KEY}: {cells_per_batch.size} batches, "
      f"min {cells_per_batch.min():,} cells ({cells_per_batch.idxmin()})", flush=True)
assert cells_per_batch.min() >= C.MIN_CELLS_PER_COHORT, (
    f"a cohort under {C.MIN_CELLS_PER_COHORT} cells reached this step "
    f"({cells_per_batch.min()}); subset_and_qc.ipynb should have dropped it"
)

# the panel
stats = hvg_stats(adata)
n_batches = adata.obs[C.BATCH_KEY].cat.categories.size
base = ranked_hvgs(stats, n_batches, C.N_HVGS)
# The ranking must match scib's, or the replacements are not "the next ones down".
scib_panel = scib.preprocessing.hvg_batch(adata, batch_key=C.BATCH_KEY,
                                           target_genes=C.N_HVGS, n_bins=20)
assert set(base) == set(scib_panel), (
    f"the reproduced ranking differs from scib's hvg_batch in "
    f"{len(set(base) ^ set(scib_panel))} genes")
print(f"Reproduced scib's {C.N_HVGS}-gene selection exactly", flush=True)

mito = [g for g in base if g.startswith(MT_PREFIX)]
target = C.N_HVGS
while True:
    # a gene just below the cut can itself be MT-, so read further down until enough remain
    panel = [g for g in ranked_hvgs(stats, n_batches, target) if not g.startswith(MT_PREFIX)]
    if len(panel) >= C.N_HVGS:
        break
    target += C.N_HVGS - len(panel)
panel = panel[: C.N_HVGS]
added = [g for g in panel if g not in set(base)]
print(f"MT- genes replaced: {len(mito)} ({', '.join(mito)})", flush=True)
print(f"replacements, the next genes of the ranking: {', '.join(added)}", flush=True)
assert len(set(panel)) == C.N_HVGS and not any(g.startswith(MT_PREFIX) for g in panel)

adata.var["highly_variable"] = adata.var_names.isin(panel)
adata.uns["hvg_panel"] = "nomt"
adata.uns["hvg_panel_dropped"] = mito
adata.uns["hvg_panel_added"] = added

# Overlap with the 01_5 HVGs, if on disk; informative only.
p = os.path.join(C.data_dir(), "shiao_hvg_2k_unintegrated_list.csv")
if os.path.exists(p):
    other = set(pd.read_csv(p, header=None)[0].astype(str))
    shared = len(set(panel) & other)
    print(f"Overlap with the 01_5 full-object HVGs: {shared}/{C.N_HVGS} "
          f"({100 * shared / C.N_HVGS:.1f}%)", flush=True)

# PCA / neighbours / UMAP
# 01_5's parameters: 50 PCs on the panel, arpack, default neighbours, UMAP.
sc.tl.pca(adata, n_comps=50, mask_var="highly_variable", svd_solver="arpack", return_info=True)
sc.pp.neighbors(adata, use_rep="X_pca")
sc.tl.umap(adata)

# `cell_type` is constant here, so the label panels use the pre-CNV label.
suffix = f"_{C.compartment()}"
label_for_plots = C.LABEL_KEY if adata.obs[C.LABEL_KEY].nunique() > 1 else C.PRIOR_LABEL_KEY
print(f"colouring the label panels by {label_for_plots!r} "
      f"({adata.obs[label_for_plots].nunique()} levels)", flush=True)

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save=f"{suffix}_elbow.png")
sc.pl.pca(adata, color=label_for_plots, save=f"{suffix}_{label_for_plots}.png")
sc.pl.umap(adata, color=label_for_plots, save=f"{suffix}_{label_for_plots}.png")
sc.pl.umap(adata, color=C.BATCH_KEY, save=f"{suffix}_{C.BATCH_KEY}.png")
if "phase" in adata.obs:
    sc.pl.umap(adata, color="phase", save=f"{suffix}_phase.png")
if "cnv_corr" in adata.obs:
    sc.pl.umap(adata, color="cnv_corr", save=f"{suffix}_cnv_corr.png")

adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)

# the DRVI input
# Built here so the .h5ad and the .var flag cannot drift apart.
print(f"\nBuilding the DRVI input ({C.N_HVGS:,} HVGs without MT-)...", flush=True)
hvg = adata[:, pd.Index(panel)].copy()
assert hvg.n_vars == C.N_HVGS, f"got {hvg.n_vars} genes after subsetting"
assert hvg.n_obs == adata.n_obs, "cells were lost during the gene subset"

# DRVI trains on .layers['counts']: fail here rather than after hours on the wrong input.
counts = hvg.layers["counts"]
counts_values = counts.data if hasattr(counts, "data") else np.asarray(counts).ravel()
assert np.allclose(counts_values, np.floor(counts_values)), (
    "layers['counts'] contains non-integer values: it is no longer raw counts"
)
assert counts_values.min() >= 0, "layers['counts'] contains negative values"

x_values = hvg.X.data if hasattr(hvg.X, "data") else np.asarray(hvg.X).ravel()
assert not np.allclose(x_values, np.floor(x_values)), (
    ".X looks like integer counts; were .X and layers['counts'] swapped upstream?"
)

hvg.write_h5ad(HVG_H5AD_PATH, compression="gzip")
print(f"Wrote {HVG_H5AD_PATH} "
      f"({hvg.n_obs:,} x {hvg.n_vars:,}, "
      f"{os.path.getsize(HVG_H5AD_PATH) / 1024 ** 3:.2f} GB on disk)", flush=True)

# The list last, so a crash never leaves a list without its object.
pd.Series(hvg.var_names).to_csv(HVG_CSV_PATH, index=False, header=False)
print(f"Wrote {len(hvg.var_names)} genes to {HVG_CSV_PATH}", flush=True)
print("\nnext: 03_3 run_drvi_tum.py trains on this panel", flush=True)
