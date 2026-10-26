#!/usr/bin/env python3
"""
03_1 prepare_infercnv_input: per-cohort inferCNV inputs from shiao.h5ad

There is no normal tissue in this dataset, so the diploid baseline comes from the immune
cells of each tumor. One run per patient, so no cell is compared against another
patient's reference.

Groups written for each cohort:
  ref_tcell       up to 1000 T/NK cells               reference
  ref_myeloid     up to 1000 myeloid cells            reference
  stromal         up to 2000 fibro/endo/PVL cells     null the per-cohort cut comes from
  immune_heldout  up to 500 immune cells              specificity control, not in the reference
  epi             every epithelial cell               cells to call
Two reference groups, so a gene higher in myeloid than in T cells is not read as a gain.
B/plasma cells are in no group (their immunoglobulin loci would add structure on chr2/14/22).

Input : $DATA_DIR/shiao.h5ad
        Read lazily from the HDF5 file. .layers['counts'] = raw counts; .obs['celltypist_predicted'] from 01_4.
Output: $DATA_DIR/03_tum/gene_order_hg38_gencode_v27.txt
        inferCNV gene ordering file (hg38, GENCODE v27), downloaded once.
        $DATA_DIR/03_tum/input/<cohort>/{counts.mtx, genes.tsv, barcodes.tsv, annotations.tsv}
        Raw counts as genes x cells Matrix Market, and the group of each cell.
        $DATA_DIR/03_tum/cohort_census.csv
        Cells per group for each cohort; infercnv_all.sh reads it to decide what to run.

Usage:
  export DATA_DIR=<data root>
  python3 prepare_infercnv_input.py
  python3 prepare_infercnv_input.py --cohorts Patient52 Patient16   # only these cohorts
  python3 prepare_infercnv_input.py --force                        # rewrite existing inputs
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from anndata.io import read_elem, sparse_dataset

# Lineage sets of 01_4/fraction_reassignment.py, split into the inferCNV groups.
# main() asserts that together they cover every cell_type label.
EPITHELIAL = {
    "LummHR-SCGB", "LummHR-active", "LummHR-major",
    "Lumsec-HLA", "Lumsec-KIT", "Lumsec-basal", "Lumsec-lac", "Lumsec-major",
    "Lumsec-myo", "Lumsec-prol",
    "basal",
}
# The stromal null is sampled stratified over these three (see stratified_sample).
STROMAL_FIBRO = {"Fibro-SFRP4", "Fibro-major", "Fibro-matrix", "Fibro-prematrix"}
STROMAL_ENDO = {
    "Lymph-immune", "Lymph-major", "Lymph-valve1", "Lymph-valve2",
    "Vas-arterial", "Vas-capillary", "Vas-venous",
}
STROMAL_PVL = {"pericytes", "vsmc"}
STROMAL_STRATA = {"fibro": STROMAL_FIBRO, "endo": STROMAL_ENDO, "pvl": STROMAL_PVL}
STROMAL = STROMAL_FIBRO | STROMAL_ENDO | STROMAL_PVL
T_NK = {
    "CD4-Tem", "CD4-Th", "CD4-Th-like", "CD4-Treg", "CD4-activated", "CD4-naive",
    "CD8-Tem", "CD8-Trm", "CD8-activated", "GD", "NK", "NK-ILCs", "NKT", "T_prol",
}
MYELOID = {
    "Macro-IFN", "Macro-lipo", "Macro-m1", "Macro-m1-CCL", "Macro-m2", "Macro-m2-CXCL",
    "Mast", "Mono-classical", "Mono-non-classical", "Neutrophil",
    "cDC1", "cDC2", "mDC", "pDC", "mye-prol",
}
# B/plasma cells are left out of every run: their immunoglobulin transcriptome would add
# fake structure on chr2/chr14/chr22.
B_PLASMA = {"b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG"}

GROUP_EPI = "epi"
GROUP_STROMAL = "stromal"
GROUP_HELDOUT = "immune_heldout"
GROUP_REF_T = "ref_tcell"
GROUP_REF_MYE = "ref_myeloid"
REFERENCE_GROUPS = (GROUP_REF_T, GROUP_REF_MYE)
# Order of the blocks in the inputs and of the census columns.
BLOCK_ORDER = (GROUP_REF_T, GROUP_REF_MYE, GROUP_STROMAL, GROUP_HELDOUT, GROUP_EPI)

BATCH_KEY = "cohort"
LABEL_KEY = "cell_type"
RAW_LABEL_KEY = "celltypist_predicted" # CellTypist before majority voting

N_REF_PER_GROUP = 1000 # T/NK and myeloid cells per cohort, each
N_STROMAL = 2000 # stromal cells in the null
N_IMMUNE_HELDOUT = 500 # immune cells kept out of the reference, as observations
MIN_REF_TOTAL = 1000 # held-out cells never leave the reference below this
MIN_EPI_CELLS = 50 # cohorts with fewer epithelial cells are not run
MIN_CELLS_PER_GENE = 3 # drop genes detected in fewer cells of the cohort
SEED = 0

GENE_ORDER_URL = "https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt"
GENE_ORDER_NAME = "gene_order_hg38_gencode_v27.txt"
MAIN_CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    p.add_argument("--cohorts", nargs="+", default=None,
                   help="only these cohorts (default: every cohort with enough epithelium)")
    p.add_argument("--force", action="store_true",
                   help="rewrite inputs for cohorts that already have them")
    return p.parse_args()


def fetch_gene_order(cnv_dir: Path) -> pd.DataFrame:
    """Download the inferCNV hg38 gene ordering file once, return it parsed."""
    path = cnv_dir / GENE_ORDER_NAME
    if not path.exists():
        print(f"Downloading gene order file -> {path}", flush=True)
        tmp = path.with_suffix(".tmp")
        urllib.request.urlretrieve(GENE_ORDER_URL, tmp)
        os.replace(tmp, path)
    go = pd.read_csv(path, sep="\t", header=None,
                     names=["gene", "chr", "start", "stop", "_trailing"],
                     usecols=[0, 1, 2, 3])
    go = go[go["chr"].isin(MAIN_CHROMS)]
    assert not go["gene"].duplicated().any(), "duplicated gene symbol in the ordering file"
    # inferCNV needs the genes sorted by position.
    go["chr"] = pd.Categorical(go["chr"], categories=MAIN_CHROMS, ordered=True)
    go = go.sort_values(["chr", "start", "stop"]).set_index("gene")
    print(f"Gene order: {len(go):,} genes on {go['chr'].nunique()} chromosomes", flush=True)
    return go


def assign_group(labels: pd.Series) -> pd.Series:
    """Map cell_type onto the inferCNV groups; B/plasma -> NaN."""
    out = pd.Series(np.nan, index=labels.index, dtype=object)
    out[labels.isin(EPITHELIAL)] = GROUP_EPI
    out[labels.isin(STROMAL)] = GROUP_STROMAL
    out[labels.isin(T_NK)] = GROUP_REF_T
    out[labels.isin(MYELOID)] = GROUP_REF_MYE
    return out


def stratified_sample(strata: dict, cap: int, rng: np.random.Generator) -> np.ndarray:
    """Draw up to cap indices, cap/len(strata) per stratum; a short stratum leaves its share to the others.

    Uniform sampling would let each cohort's stromal mix define its null, and the three
    populations do not score at the same level.
    """
    share = cap // len(strata)
    taken = {}
    for name, idx in strata.items():
        n = min(share, len(idx))
        taken[name] = rng.choice(idx, size=n, replace=False) if n else idx[:0]
    short = cap - sum(len(v) for v in taken.values())
    if short > 0:
        left = np.concatenate([np.setdiff1d(idx, taken[name])
                               for name, idx in strata.items()])
        if len(left):
            taken["_topup"] = rng.choice(left, size=min(short, len(left)), replace=False)
    parts = [v for v in taken.values() if len(v)]
    if not parts:
        return np.array([], dtype=int)
    return np.sort(np.concatenate(parts).astype(int))


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(SEED)

    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    in_path = data_dir / "shiao.h5ad"
    cnv_dir = data_dir / "03_tum"
    input_dir = cnv_dir / "input"
    cnv_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(exist_ok=True)

    print(f"DATA_DIR: {data_dir}")
    print(f"Input   : {in_path} ({in_path.stat().st_size / 1024**3:.2f} GB)")
    print(f"Output  : {input_dir}")
    print()

    gene_order = fetch_gene_order(cnv_dir)

    # Here only .obs/.var are read, and one cohort's counts at a time.
    h5 = h5py.File(in_path, "r")
    assert "counts" in h5["layers"], "expected raw counts in layers['counts']"
    obs = read_elem(h5["obs"])
    var_names = read_elem(h5["var"]).index
    counts_ds = sparse_dataset(h5["layers"]["counts"])

    assert RAW_LABEL_KEY in obs, (
        f"obs['{RAW_LABEL_KEY}'] is needed to keep tumor cells out of the reference; "
        "it is written by 01_4/celltypist_annotation.py"
    )
    labels = obs[LABEL_KEY].astype(str)
    raw_labels = obs[RAW_LABEL_KEY].astype(str)
    unclassified = set(labels.unique()) - EPITHELIAL - STROMAL - T_NK - MYELOID - B_PLASMA
    assert not unclassified, (
        f"cell_type value(s) in no lineage set: {sorted(unclassified)}; "
        "update the sets at the top of this script before running"
    )
    groups = assign_group(labels)

    # Genes intersected once on the full object, so every cohort starts from the same genes.
    keep_genes = var_names.intersection(gene_order.index)
    keep_genes = pd.Index(gene_order.index[gene_order.index.isin(keep_genes)]) # positional order
    gene_pos = pd.Series(np.arange(len(var_names)), index=var_names)[keep_genes].to_numpy()
    print(f"Genes: {len(var_names):,} in the object, {len(keep_genes):,} also in the gene order file")
    print()

    cohorts = list(obs[BATCH_KEY].cat.categories) if args.cohorts is None else args.cohorts
    census_rows = []

    for cohort in cohorts:
        in_cohort = (obs[BATCH_KEY].astype(str) == cohort).to_numpy()
        grp = groups[in_cohort]
        n_epi = int((grp == GROUP_EPI).sum())

        if n_epi < MIN_EPI_CELLS:
            print(f"[skip] {cohort}: {n_epi} epithelial cells (< {MIN_EPI_CELLS})")
            census_rows.append(dict(cohort=cohort, n_epi=n_epi, status="skipped_too_few_epi"))
            continue

        out_dir = input_dir / cohort
        annot_path = out_dir / "annotations.tsv"
        if annot_path.exists() and not args.force:
            # Still add it to the census, which infercnv_all.sh reads to decide what to run.
            have = pd.read_csv(annot_path, sep="\t", header=None, names=["cell", "group"])
            have_n = have["group"].value_counts()
            print(f"[have] {cohort}: inputs already written, skipping")
            census_rows.append(dict(
                cohort=cohort, status="prepared",
                n_cells=len(have),
                n_genes=sum(1 for _ in open(out_dir / "genes.tsv")),
                **{f"n_{g}": int(have_n.get(g, 0)) for g in BLOCK_ORDER},
            ))
            continue

        # Held-out immune cells are drawn before the reference and removed from its pool.
        labels_c = labels[in_cohort]
        raw_c = raw_labels[in_cohort]
        picked = {GROUP_EPI: np.flatnonzero((grp == GROUP_EPI).to_numpy())}

        strata = {name: np.flatnonzero(labels_c.isin(members).to_numpy())
                  for name, members in STROMAL_STRATA.items()}
        picked[GROUP_STROMAL] = stratified_sample(strata, N_STROMAL, rng)

        # Immune cells whose raw (pre-voting) CellTypist call is epithelial are excluded:
        # majority voting can give a tumor cell an immune label.
        pure = ~raw_c.isin(EPITHELIAL).to_numpy()
        is_imm = ((grp == GROUP_REF_T) | (grp == GROUP_REF_MYE)).to_numpy()
        imm = np.flatnonzero(is_imm & pure)
        n_held = min(N_IMMUNE_HELDOUT, max(0, len(imm) - MIN_REF_TOTAL))
        held = np.sort(rng.choice(imm, size=n_held, replace=False)) if n_held else imm[:0]
        picked[GROUP_HELDOUT] = held

        avail = np.setdiff1d(imm, held)
        for group in REFERENCE_GROUPS:
            idx = np.intersect1d(avail, np.flatnonzero((grp == group).to_numpy()))
            if len(idx) > N_REF_PER_GROUP:
                idx = rng.choice(idx, size=N_REF_PER_GROUP, replace=False)
            picked[group] = np.sort(idx)

        # Positions are relative to the cohort slice; map them back to the full object.
        cohort_pos = np.flatnonzero(in_cohort)
        order = np.concatenate([picked[g] for g in BLOCK_ORDER])
        rows = cohort_pos[order]
        cell_groups = np.concatenate([[g] * len(picked[g]) for g in BLOCK_ORDER])
        cell_ids = obs.index[rows]
        strat_n = {f"n_stromal_{k}": int(np.isin(picked[GROUP_STROMAL], v).sum())
                   for k, v in strata.items()}

        # Backed slicing needs increasing indices: sort the rows and reorder the labels to match.
        srt = np.argsort(rows, kind="stable")
        counts = counts_ds[rows[srt]]
        counts = sp.csr_matrix(counts)[:, gene_pos]
        cell_ids = cell_ids[srt]
        cell_groups = cell_groups[srt]

        # Genes barely detected in this cohort carry no CNV signal.
        detected = np.asarray((counts > 0).sum(axis=0)).ravel()
        gene_mask = detected >= MIN_CELLS_PER_GENE
        counts = counts[:, gene_mask]
        genes_here = keep_genes[gene_mask]

        assert counts.shape == (len(cell_ids), len(genes_here))
        n_by_group = pd.Series(cell_groups).value_counts()

        out_dir.mkdir(parents=True, exist_ok=True)
        # inferCNV wants genes as rows and cells as columns.
        sio.mmwrite(str(out_dir / "counts.mtx"), counts.T.astype(np.int32), field="integer")
        pd.Series(genes_here).to_csv(out_dir / "genes.tsv", index=False, header=False)
        pd.Series(cell_ids).to_csv(out_dir / "barcodes.tsv", index=False, header=False)
        # inferCNV annotation format: cell id <tab> group, no header.
        pd.DataFrame({"cell": cell_ids, "group": cell_groups}).to_csv(
            out_dir / "annotations.tsv", sep="\t", index=False, header=False)

        mtx_mb = (out_dir / "counts.mtx").stat().st_size / 1024**2
        print(f"[ok]   {cohort}: {counts.shape[0]:,} cells x {counts.shape[1]:,} genes "
              f"({mtx_mb:.0f} MB) | " +
              " ".join(f"{g}={n_by_group.get(g, 0)}" for g in BLOCK_ORDER) +
              " | null: " + " ".join(f"{k.split('_')[-1]}={v}" for k, v in strat_n.items()),
              flush=True)

        census_rows.append(dict(
            cohort=cohort, status="prepared",
            n_cells=counts.shape[0], n_genes=counts.shape[1],
            **{f"n_{g}": int(n_by_group.get(g, 0)) for g in BLOCK_ORDER},
            # composition of the stromal null
            **strat_n,
        ))

    census = pd.DataFrame(census_rows)
    census_path = cnv_dir / "cohort_census.csv"
    if census_path.exists() and args.cohorts is not None:
        # A partial re-run keeps the census rows of the other cohorts.
        old = pd.read_csv(census_path)
        census = pd.concat([old[~old["cohort"].isin(census["cohort"])], census])
    census.sort_values("cohort").to_csv(census_path, index=False)
    print()
    print(f"Wrote {census_path}")
    prepared = census[census["status"] == "prepared"] if len(census) else census
    print(f"{len(prepared)} cohort(s) ready for run_infercnv.R")
    return 0


if __name__ == "__main__":
    sys.exit(main())
