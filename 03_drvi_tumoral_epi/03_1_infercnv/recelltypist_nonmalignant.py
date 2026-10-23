#!/usr/bin/env python3
"""
03_1 recelltypist_nonmalignant: re-run the 01_4 CellTypist annotation without the malignant cells

Malignant cells get the single label 'malignant' and every other cell
 is re-annotated with the same model and settings as 01_4.

All non-malignant cells are re-annotated: the epithelial set is defined by the labels
being repaired, and a smaller input would change the voting neighbourhoods. 
Labels also move through the stochastic over-clustering (e.g. CD4 <-> CD8),
so the label churn is not a measure of contamination.

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

import celltypist
from celltypist import models
from anndata.io import read_elem, sparse_dataset

sc.settings.verbosity = 1

MALIGNANT_LABEL = "malignant"
STATUS_KEY = "cnv_status"
LABEL_KEY = "cell_type"
BATCH_KEY = "cohort"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    p.add_argument("--no-majority-voting", action="store_true",
                   help="skip CellTypist's over-clustering vote (fast, but not what 01_4 did)")
    p.add_argument("--scope", choices=("all", "non-immune"), default="all",
                   help="non-malignant cells to re-annotate; 'all' (default) or only the "
                        "non-immune fraction (faster)")
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

    print("Building temporary CP10K+log1p matrix for CellTypist input...", flush=True)
    # No phase-01 PCA/graph/UMAP is carried over, so the over-clustering is rebuilt without the
    # tumour cells.
    rows = np.flatnonzero(reannotate)
    counts = sparse_dataset(h5["layers"]["counts"])[rows]
    h5.close()
    adata_ct = ad.AnnData(X=counts, obs=obs.iloc[rows].copy(), var=var)
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    majority_voting = not args.no_majority_voting
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
