"""
01_5 scib reduce_data
Feature selection (batch-aware HVG) + PCA + neighbors graph + UMAP on the
full, already scran-normalized dataset. No (re)normalization happens here:
.X is consumed as-is from the 01_4 output.

Input : $DATA_DIR/all_samples_combined_scrublet_norm_cc_annotated.h5ad
        .X       = scran log-normalized
        .layers['counts'] = raw integer counts
        .var     = gene symbols as var_names, Ensembl in .var['gene_ids']
        .obs     = 'cohort' (34 patients, integration batch_key), 'cell_type'
                   (CellTypist majority_voting label_key), etc.
Output: $DATA_DIR/all_samples_combined_scrublet_norm_cc_annotated_reduced.h5ad
        adds .var['highly_variable'] (binary, batch-aware), .obsm['X_pca'],
        neighbors graph, .obsm['X_umap']
        $DATA_DIR/shiao_hvg_2k_unintegrated_list.csv
        gene symbol + Ensembl id of the selected HVGs (one row per gene)

Usage:
  export DATA_DIR=<data root>
  python3 scib_reduce_data.py
"""

import os
import matplotlib
matplotlib.use("Agg")  # headless backend: figures are only written to disk via save=, never shown/blocking
import numpy as np
import pandas as pd

# Fix for an error in scib's hvg_batch: it calls np.in1d, removed in numpy >= 2.4.
# Map it to np.isin.
if not hasattr(np, "in1d"):
    np.in1d = np.isin

import anndata2ri
from rpy2.robjects import conversion, default_converter

def _activate():
    conversion.set_conversion(default_converter + anndata2ri.converter)

def _deactivate():
    conversion.set_conversion(default_converter)

anndata2ri.activate = _activate
anndata2ri.deactivate = _deactivate

import scanpy as sc
import scib

# Resolve DATA_DIR from environment, never hardcoded
DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm_cc_annotated.h5ad")
OUT_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm_cc_annotated_reduced.h5ad")

# Figure dir anchored to the script location, fallback to cwd for interactive use
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()
PHASE_DIR = os.path.dirname(SCRIPT_DIR)  # 01_pre_processing/
FIG_DIR = os.path.join(PHASE_DIR, "figures", "01_5_scib_reduce_data")
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=300, facecolor="white")

adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

# 2000 HVGs selected per patient (cohort) and merged, scib default settings;
# then PCA (50 comps), neighbors and UMAP.
scib.preprocessing.reduce_data(
    adata,
    batch_key="cohort",
    flavor="cell_ranger",
    n_top_genes=2000,
    n_bins=20,
    pca=True,
    pca_comps=50,
    svd_solver="arpack",
    overwrite_hvg=True,
    neighbors=True,
    use_rep="X_pca",
    umap=True,
)
print(f"HVGs selected: {int(adata.var['highly_variable'].sum())}")

# Save the selected HVGs for reuse outside this AnnData (e.g. restricting other objects/tools
# to the same feature set): plain list of gene symbols, one per line, no header.
hvg_genes = adata.var_names[adata.var["highly_variable"]]
HVG_CSV_PATH = os.path.join(DATA_DIR, "shiao_hvg_2k_unintegrated_list.csv")
pd.Series(hvg_genes).to_csv(HVG_CSV_PATH, index=False, header=False)
print(f"Wrote {len(hvg_genes)} HVGs to {HVG_CSV_PATH}")

# Diagnostic figures
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save="_elbow.png")
sc.pl.pca(adata, color="fraction", save="_fraction.png")
sc.pl.pca(adata, color="cell_type", save="_cell_type.png")
sc.pl.umap(adata, color="cell_type", save="_cell_type.png") 

adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}")