#!/usr/bin/env python3
"""
03_1 recelltypist_nonmalignant: re-run the 01_4 CellTypist annotation without the malignant cells

Malignant cells get the single label 'malignant' and every other cell
 is re-annotated with the same model and settings as 01_4.

All non-malignant cells are re-annotated: the epithelial set is defined by the labels
being repaired, and a smaller input would change the voting neighbourhoods. 

--low-memory runs the same annotation in blocks of cells (~7 GB instead of >30 GB).

Input : $DATA_DIR/shiao.h5ad
        .layers['counts'] = raw counts.
        $DATA_DIR/03_tum/cnv_status.csv
        Output of call_malignant.ipynb.
        $DATA_DIR/Cells_Adult_Breast.pkl
        Same model as 01_4.
Output: $DATA_DIR/03_tum/cell_annotation.csv
        One row per cell: cohort, cell_type (01_4), cnv_status, cell_type_cnv (new label,
        'malignant' for malignant cells), celltypist_predicted_cnv (before voting, for QC),
        reannotated. borderline and not_tested cells also get a normal-breast label, so
        03_2 must read cnv_status, not the label.

Usage:
  export DATA_DIR=<data root>
  python3 recelltypist_nonmalignant.py
  python3 recelltypist_nonmalignant.py --no-majority-voting    # per-cell only, much faster
  python3 recelltypist_nonmalignant.py --scope non-immune      # only fraction == 'non_imm'
  python3 recelltypist_nonmalignant.py --low-memory            # in blocks, fits in ~7 GB
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

import celltypist
from celltypist import models
from celltypist.classifier import AnnotationResult, Classifier
from anndata.io import read_elem, sparse_dataset
# scanpy's own seurat-HVG steps, so --low-memory picks the genes exactly as sc.pp.highly_variable_genes
from scanpy.preprocessing._highly_variable_genes import _get_disp_stats, _get_mean_bins, _subset_genes

sc.settings.verbosity = 1

MALIGNANT_LABEL = "malignant"
STATUS_KEY = "cnv_status"
LABEL_KEY = "cell_type"
BATCH_KEY = "cohort"

# Classifier._construct_neighbor_graph settings
GRAPH_MIN_CELLS = 5
GRAPH_N_HVG = 2500
GRAPH_MAX_VALUE = 10
GRAPH_N_PCS = 50
GRAPH_N_NEIGHBORS = 10


def lognorm(counts_ds, rows: np.ndarray, obs: pd.DataFrame, var: pd.DataFrame) -> ad.AnnData:
    """CP10K + log1p of the raw counts of rows, the input CellTypist expects."""
    a = ad.AnnData(X=counts_ds[rows], obs=obs.iloc[rows].copy(), var=var)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def annotate_low_memory(counts_ds, rows: np.ndarray, obs: pd.DataFrame, var: pd.DataFrame,
                        model, majority_voting: bool, block_size: int) -> tuple[pd.Series, pd.Series]:
    """celltypist.annotate in blocks of cells; returns (per-cell labels, voted labels).

    The over-clustering repeats Classifier._construct_neighbor_graph from per-gene sums and a
    covariance accumulated over blocks, then Classifier.over_cluster and majority_vote.
    """
    blocks = [rows[i:i + block_size] for i in range(0, len(rows), block_size)]
    n, n_genes = len(rows), var.shape[0]

    # Pass 1: per-cell prediction, and the per-gene sums for the gene filter and seurat HVGs
    raw = []
    ncell, s1, s2 = np.zeros(n_genes), np.zeros(n_genes), np.zeros(n_genes)
    for i, b in enumerate(blocks, 1):
        a = lognorm(counts_ds, b, obs, var)
        pred = celltypist.annotate(a, model=model, majority_voting=False)
        raw.append(pred.predicted_labels["predicted_labels"].astype(str))
        if majority_voting:
            x = a.X.astype(np.float64)
            ncell += np.asarray((x > 0).sum(axis=0)).ravel()
            x.data = np.expm1(x.data)  # seurat dispersion is on the linear scale
            s1 += np.asarray(x.sum(axis=0)).ravel()
            s2 += np.asarray(x.multiply(x).sum(axis=0)).ravel()
        print(f"  block {i}/{len(blocks)}: predicted", flush=True)
    raw = pd.concat(raw)
    if not majority_voting:
        return raw, raw

    # Genes in >= 5 cells, then the 2,500 seurat HVGs, as in sc.pp.highly_variable_genes
    keep = ncell >= GRAPH_MIN_CELLS
    mean = s1[keep] / n
    disp = (s2[keep] - n * mean**2) / (n - 1)
    mean[mean == 0] = 1e-12
    disp /= mean
    disp[disp == 0] = np.nan
    df = pd.DataFrame({"means": np.log1p(mean), "dispersions": np.log(disp)})
    df["mean_bin"] = _get_mean_bins(df["means"], "seurat", 20)
    stats = _get_disp_stats(df, "seurat")
    disp_norm = ((df["dispersions"] - stats["avg"]) / stats["dev"]).to_numpy()
    is_hvg = _subset_genes(ad.AnnData(np.zeros((1, int(keep.sum())))), mean=df["means"].to_numpy(),
                           dispersion_norm=disp_norm, cutoff=min(GRAPH_N_HVG, int(keep.sum())))
    hvg = np.flatnonzero(keep)[is_hvg]
    print(f"  {len(hvg):,} HVGs for the over-clustering graph", flush=True)

    # Pass 2: the HVG columns of every cell, kept sparse (a fraction of the full matrix)
    h = sp.vstack([lognorm(counts_ds, b, obs, var).X[:, hvg] for b in blocks]).tocsr()
    slices = [slice(i, min(i + block_size, n)) for i in range(0, n, block_size)]

    def dense(sl):
        return h[sl].toarray().astype(np.float64)

    # sc.pp.scale(max_value=10): zero-centred, unit variance (ddof=1), clipped
    m1, m2 = np.zeros(len(hvg)), np.zeros(len(hvg))
    for sl in slices:
        z = dense(sl); m1 += z.sum(axis=0); m2 += (z**2).sum(axis=0)
    mu = m1 / n
    sd = np.sqrt((m2 - n * mu**2) / (n - 1))
    sd[sd == 0] = 1

    def scaled(sl):
        z = dense(sl); z -= mu; z /= sd
        return np.clip(z, -GRAPH_MAX_VALUE, GRAPH_MAX_VALUE, out=z)

    # sc.pp.pca(n_comps=50) from the covariance: same components as the dense PCA, up to sign
    z1, ztz = np.zeros(len(hvg)), np.zeros((len(hvg), len(hvg)))
    for sl in slices:
        z = scaled(sl); z1 += z.sum(axis=0); ztz += z.T @ z
    zm = z1 / n
    _, vecs = np.linalg.eigh((ztz - n * np.outer(zm, zm)) / (n - 1))
    vecs = vecs[:, ::-1][:, :GRAPH_N_PCS]
    pcs = np.vstack([(scaled(sl) - zm) @ vecs for sl in slices]).astype(np.float32)
    del h

    g = ad.AnnData(obs=pd.DataFrame(index=raw.index))
    g.obsm["X_pca"] = pcs
    sc.pp.neighbors(g, n_neighbors=GRAPH_N_NEIGHBORS, n_pcs=GRAPH_N_PCS, use_rep="X_pca")
    # over_cluster only needs .adata: it finds the graph and picks the resolution by cell count
    clf = Classifier.__new__(Classifier)
    clf.adata = g
    over = clf.over_cluster()
    voted = Classifier.majority_vote(
        AnnotationResult(raw.to_frame("predicted_labels"), None, None, None), over)
    return raw, voted.predicted_labels["majority_voting"].astype(str)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    p.add_argument("--no-majority-voting", action="store_true",
                   help="skip CellTypist's over-clustering vote (fast, but not what 01_4 did)")
    p.add_argument("--scope", choices=("all", "non-immune"), default="all",
                   help="non-malignant cells to re-annotate; 'all' (default) or only the "
                        "non-immune fraction (faster)")
    p.add_argument("--low-memory", action="store_true",
                   help="annotate in blocks of cells (~7 GB instead of >30 GB); per-cell labels "
                        "identical, voted labels equivalent but not bit-identical")
    p.add_argument("--block-size", type=int, default=50_000,
                   help="cells per block with --low-memory [%(default)s]")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    in_path = data_dir / "shiao.h5ad"
    cnv_dir = data_dir / "03_tum"
    status_path = cnv_dir / "cnv_status.csv"
    model_path = data_dir / "Cells_Adult_Breast.pkl"
    out_path = cnv_dir / "cell_annotation.csv"

    for path, hint in ((in_path, "phase 01"),
                       (status_path, "call_malignant.ipynb"),
                       (model_path, "01_4")):
        assert path.exists(), f"missing {path} (produced by {hint})"

    print(f"DATA_DIR: {data_dir}")
    print(f"Input   : {in_path}")
    print(f"Status  : {status_path}")
    print(f"Output  : {out_path}")
    print()

    status = pd.read_csv(status_path, index_col=0)
    assert STATUS_KEY in status.columns, f"{status_path} has no '{STATUS_KEY}' column"

    # Only .obs/.var and the counts of the cells to re-annotate are read: the whole object
    # (.X + counts, ~15 GB in memory) does not fit on a workstation.
    print("Loading .obs...", flush=True)
    h5 = h5py.File(in_path, "r")
    assert "counts" in h5["layers"], "expected raw counts in layers['counts']"
    obs = read_elem(h5["obs"])
    var = read_elem(h5["var"])

    # The call has one row per cell of shiao.h5ad; a missing row means the two are out of sync.
    missing = obs.index.difference(status.index)
    assert len(missing) == 0, (
        f"{len(missing):,} cells of shiao.h5ad have no row in {status_path.name}; "
        "re-run call_malignant.ipynb against the same object"
    )
    obs[STATUS_KEY] = status.loc[obs.index, STATUS_KEY].to_numpy()
    print(obs[STATUS_KEY].value_counts().to_string())
    print()

    is_malignant = (obs[STATUS_KEY] == "malignant").to_numpy()

    # Cells left out keep their 01_4 label (reannotated = False in the output).
    if args.scope == "all":
        reannotate = ~is_malignant
    else:
        assert "fraction" in obs, "obs['fraction'] is needed by --scope non-immune"
        reannotate = (~is_malignant) & (obs["fraction"] == "non_imm").to_numpy()

    print(f"malignant: {is_malignant.sum():,} cells -> single label '{MALIGNANT_LABEL}'")
    print(f"re-annotated (--scope {args.scope}): {reannotate.sum():,} cells")
    kept_as_is = (~is_malignant) & (~reannotate)
    if kept_as_is.any():
        print(f"keeping the 01_4 label unchanged: {kept_as_is.sum():,} cells")
    print()

    # Same annotation as 01_4: temporary CP10K+log1p copy from raw counts.
    print(f"Loading CellTypist model from: {model_path}", flush=True)
    model = models.Model.load(model=str(model_path))

    # No phase-01 PCA/graph/UMAP is carried over, so the over-clustering is rebuilt without the
    # tumour cells. 01_5's HVG column would replace CellTypist's own gene choice (see top).
    rows = np.flatnonzero(reannotate)
    var = var.drop(columns=["highly_variable"], errors="ignore")
    majority_voting = not args.no_majority_voting
    counts_ds = sparse_dataset(h5["layers"]["counts"])

    if args.low_memory:
        print(f"Annotating with CellTypist in blocks of {args.block_size:,} cells "
              f"(majority_voting={majority_voting}, CPU)...", flush=True)
        raw_labels, new_labels = annotate_low_memory(counts_ds, rows, obs, var, model,
                                                     majority_voting, args.block_size)
        h5.close()
    else:
        print("Building temporary CP10K+log1p matrix for CellTypist input...", flush=True)
        adata_ct = lognorm(counts_ds, rows, obs, var)
        h5.close()

        print(f"Annotating with CellTypist (majority_voting={majority_voting}, CPU)...", flush=True)
        predictions = celltypist.annotate(
            adata_ct,
            model=model,
            majority_voting=majority_voting,
            use_GPU=False,
        )
        voted_col = "majority_voting" if majority_voting else "predicted_labels"
        new_labels = predictions.predicted_labels[voted_col].astype(str)
        raw_labels = predictions.predicted_labels["predicted_labels"].astype(str)
        del adata_ct

    # Per-cell table
    out = pd.DataFrame(index=obs.index)
    out[BATCH_KEY] = obs[BATCH_KEY].astype(str).to_numpy()
    out[LABEL_KEY] = obs[LABEL_KEY].astype(str).to_numpy()
    out[STATUS_KEY] = obs[STATUS_KEY].to_numpy()

    # 01_4 label, then 'malignant' for malignant cells, then the new labels for re-annotated cells.
    out["cell_type_cnv"] = out[LABEL_KEY]
    out["celltypist_predicted_cnv"] = out[LABEL_KEY]
    out.loc[is_malignant, "cell_type_cnv"] = MALIGNANT_LABEL
    out.loc[is_malignant, "celltypist_predicted_cnv"] = MALIGNANT_LABEL
    out.loc[new_labels.index, "cell_type_cnv"] = new_labels.to_numpy()
    out.loc[raw_labels.index, "celltypist_predicted_cnv"] = raw_labels.to_numpy()
    out["reannotated"] = reannotate

    assert (out.loc[is_malignant, "cell_type_cnv"] == MALIGNANT_LABEL).all()
    assert (out.loc[~is_malignant, "cell_type_cnv"] != MALIGNANT_LABEL).all(), \
        "the model produced a label literally called 'malignant'; rename MALIGNANT_LABEL"
    assert out["cell_type_cnv"].notna().all()

    print()
    print("cell_type_cnv distribution:")
    print(out["cell_type_cnv"].value_counts().to_string())
    print()
    print("How the re-annotation moved the non-malignant cells (01_4 -> now), top 20 changes:")
    nm = out.loc[reannotate]
    changed = nm[nm[LABEL_KEY] != nm["cell_type_cnv"]]
    print(f"{len(changed):,} of {len(nm):,} re-annotated cells changed label "
          f"({len(changed) / max(len(nm), 1):.1%})")
    print(changed.groupby([LABEL_KEY, "cell_type_cnv"]).size()
          .sort_values(ascending=False).head(20).to_string())

    out.index.name = "cell"
    out.to_csv(out_path)
    print()
    print(f"Wrote {out_path} ({len(out):,} cells)")
    print("next: 03_2_subsetting reads this table instead of `cell_type` from shiao.h5ad")
    return 0


if __name__ == "__main__":
    sys.exit(main())
