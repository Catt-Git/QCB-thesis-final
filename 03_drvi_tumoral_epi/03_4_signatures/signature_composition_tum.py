#!/usr/bin/env python3
"""03_4b: what is actually inside each signature score - which genes carry it, and are they measurable.

A list can map most of its symbols and still score as a proxy for one gene: score_genes
averages the mapped genes, so a highly expressed gene dominates and undetected ones add zeros.
Per gene: detection_rate, share_of_mean_signal, share_of_score_variance (Cov(x_g / n, m) / Var(m))
and the Spearman with the 03_6 score. Per signature: effective_n_genes (inverse Simpson) and
frac_cells_all_zero, plus three flags (dominated_by_one_gene, low_effective_n,
undetectable_in_many_cells). Diagnostic only: nothing is dropped.

Runs after 03_6 to join its scores; without them the Spearman column is left empty.

Usage:
    export DATA_DIR=<data root>
    python signature_composition_tum.py                      # the scie collection, the default
    python signature_composition_tum.py --collection emt     # the same, on the EMT lists
    python signature_composition_tum.py --collection gavish  # the 22 TNBC-relevant metaprograms
    python signature_composition_tum.py --collection gavish --all-metaprograms  # all 41
"""

from __future__ import annotations

import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402


# Flag thresholds, not from the literature. One gene above this share of the score variance
# makes the rest of the list decoration.
DOMINANT_GENE_SHARE = 0.30

# Same floor as the mapped count in 03_4.
MIN_EFFECTIVE_N_GENES = C.MIN_SIGNATURE_GENES

# effective_n_fraction is reported, not flagged: every long list sits at 0.18-0.28 anyway.

# Above this share of all-zero cells the score mostly reports its control set.
MAX_ALL_ZERO_FRACTION = 0.10

# Genes under this detection rate are at the droplet noise floor. Counted, not dropped.
NOISE_FLOOR_DETECTION = 0.01


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    p.add_argument("--no-spearman", action="store_true",
                   help="skip the per-gene correlation with the 03_6 score (the slow part)")
    return p.parse_args()


def spearman_against(scores: np.ndarray, X: sp.spmatrix, chunk: int = 250) -> np.ndarray:
    """Spearman of every column of `X` against `scores`, in column chunks (the dense ranks are large)."""
    y = rankdata(scores)
    y = (y - y.mean()) / y.std()
    out = np.empty(X.shape[1], dtype=float)
    for start in range(0, X.shape[1], chunk):
        block = X[:, start:start + chunk].toarray()
        r = np.apply_along_axis(rankdata, 0, block)
        sd = r.std(axis=0)
        r = (r - r.mean(axis=0)) / np.where(sd > 0, sd, 1.0)
        out[start:start + chunk] = np.where(sd > 0, (y[:, None] * r).mean(axis=0), np.nan)
    return out


def gini(shares: np.ndarray) -> float:
    """Gini of the (non-negative) contribution shares. 0 = every gene contributes equally."""
    v = np.sort(np.asarray(shares, dtype=float))
    n = len(v)
    if n == 0 or v.sum() <= 0:
        return np.nan
    return float((2 * np.arange(1, n + 1) - n - 1).dot(v) / (n * v.sum()))


def main():
    args = parse_args()
    coll = SC.resolve(args)

    C.banner(f"03_4b - signature composition: {coll.title}")
    print(f"question    which genes carry each score, and are they measurable in this object?")
    print(f"object      {C.FULL_H5AD}")

    # The .gmt holds the mapped genes, the lists as Route A scores them.
    mapped = C.read_gmt(C.gmt_path(coll))
    names = [n for n in coll.names if n in mapped]
    union = sorted({g for n in names for g in mapped[n]})
    print(f"signatures  {len(names)} lists, {len(union):,} distinct mapped genes")

    # ------------------------------------------------------------------ scores
    # Joined when 03_6 has run; its absence is reported, not fatal.
    scores = None
    if args.no_spearman:
        print("[note] --no-spearman: skipping the per-gene correlation with the 03_6 score")
    elif C.scores_csv(coll).exists():
        scores = pd.read_csv(C.scores_csv(coll), index_col=0, comment="#")
        print(f"scores      {C.scores_csv(coll).name}")
    else:
        print(f"[note] {C.scores_csv(coll).name} not found - run 03_6 first for the "
              f"spearman_with_score column. Everything else is computed without it.")

    # ------------------------------------------------------------------ matrix
    adata = ad.read_h5ad(C.FULL_H5AD)
    missing = [g for g in union if g not in adata.var_names]
    if missing:                       # cannot happen from a .gmt this phase wrote; assert anyway
        sys.exit(f"[STOP] {len(missing)} genes of the .gmt are not in {C.FULL_H5AD.name}: "
                 f"{', '.join(missing[:5])} ... - the .gmt and the object have drifted apart.")
    X = adata[:, union].X
    X = sp.csc_matrix(X) if not sp.issparse(X) else X.tocsc()
    n_cells = X.shape[0]
    col_of = {g: i for i, g in enumerate(union)}
    if scores is not None:
        if list(scores.index) != list(adata.obs_names):
            sys.exit("[STOP] the 03_6 score table is not in the object's cell order.")
    print(f"cells       {n_cells:,} x {len(union):,} genes (log-normalised, all-genes object)")

    # ------------------------------------------------------------ per gene rows
    per_gene, summary = [], []
    for name in names:
        genes = mapped[name]
        idx = [col_of[g] for g in genes]
        Xs = X[:, idx].tocsc()
        n = len(genes)

        detection = np.asarray((Xs > 0).sum(axis=0)).ravel() / n_cells
        mean_expr = np.asarray(Xs.mean(axis=0)).ravel()

        total = mean_expr.sum()
        share_mean = mean_expr / total if total > 0 else np.full(n, np.nan)

        # Cov(x_g / n, m) / Var(m) sums to 1 over the list; X^T m keeps it sparse.
        m = np.asarray(Xs.mean(axis=1)).ravel()
        var_m = m.var()
        if var_m > 0:
            cov = (Xs.T @ (m - m.mean())) / n_cells
            cov = np.asarray(cov).ravel()
            share_var = cov / (n * var_m)
        else:
            share_var = np.full(n, np.nan)

        all_zero = float((np.asarray((Xs > 0).sum(axis=1)).ravel() == 0).mean())

        rho = np.full(n, np.nan)
        if scores is not None and f"score_{name}" in scores.columns:
            rho = spearman_against(scores[f"score_{name}"].to_numpy(), Xs)

        order = np.argsort(-np.nan_to_num(share_var, nan=-np.inf))
        rank_of = np.empty(n, dtype=int)
        rank_of[order] = np.arange(1, n + 1)

        per_gene.append(pd.DataFrame({
            "signature": name,
            "axis": coll.axis_of[name],
            "gene": genes,
            "detection_rate": detection,
            "mean_expression": mean_expr,
            "share_of_mean_signal": share_mean,
            "share_of_score_variance": share_var,
            "spearman_with_score": rho,
            "rank_by_variance_share": rank_of,
        }))

        eff_n = float(1.0 / np.square(share_mean).sum()) if total > 0 else np.nan
        top = order[0]
        top3_var = float(np.nan_to_num(share_var, nan=0.0)[order[:3]].sum())
        flags = []
        if share_var[top] > DOMINANT_GENE_SHARE:
            flags.append("dominated_by_one_gene")
        if eff_n < MIN_EFFECTIVE_N_GENES:
            flags.append("low_effective_n")
        if all_zero > MAX_ALL_ZERO_FRACTION:
            flags.append("undetectable_in_many_cells")

        summary.append({
            "signature": name,
            "axis": coll.axis_of[name],
            "primary": name in coll.primary_names(),
            "n_mapped": n,
            "effective_n_genes": eff_n,
            "effective_n_fraction": eff_n / n if n else np.nan,
            "top_gene": genes[top],
            "top_gene_share_of_variance": float(share_var[top]),
            "top_gene_share_of_mean_signal": float(share_mean[top]),
            "top_gene_detection_rate": float(detection[top]),
            "top_gene_spearman_with_score": float(rho[top]),
            "share_top3_of_variance": top3_var,
            "median_detection_rate": float(np.median(detection)),
            "n_genes_below_noise_floor": int((detection < NOISE_FLOOR_DETECTION).sum()),
            "frac_genes_below_noise_floor": float((detection < NOISE_FLOOR_DETECTION).mean()),
            "frac_cells_all_zero": all_zero,
            "gini_of_mean_signal": gini(share_mean),
            "flags": ",".join(flags) if flags else "none",
        })

    genes_tbl = pd.concat(per_gene, ignore_index=True)
    genes_tbl = genes_tbl.sort_values(["signature", "rank_by_variance_share"])
    summary_tbl = pd.DataFrame(summary).set_index("signature")
    summary_tbl = summary_tbl.loc[[n for n in coll.order(names) if n in summary_tbl.index]]

    # ----------------------------------------------------------------- reports
    print("\nconcentration: how many genes each score behaves as if it had")
    print(summary_tbl[["axis", "n_mapped", "effective_n_genes", "effective_n_fraction",
                       "top_gene", "top_gene_share_of_variance", "frac_cells_all_zero",
                       "flags"]].to_string(float_format="%.3f"))

    print("\nthe three genes that carry each score (share of its variance across cells)")
    for name in summary_tbl.index:
        top = genes_tbl[genes_tbl["signature"] == name].head(3)
        parts = [f"{r.gene} {r.share_of_score_variance:+.1%} (det {r.detection_rate:.1%})"
                 for r in top.itertuples()]
        print(f"  {name:24s} {'  |  '.join(parts)}")

    flagged = summary_tbl.index[summary_tbl["flags"] != "none"].tolist()
    if flagged:
        print(f"\n[warn] {len(flagged)} list(s) flagged - their score should not be read as the "
              f"biology the name claims without saying so:")
        for name in flagged:
            print(f"  {name:24s} {summary_tbl.loc[name, 'flags']}")
        print("  Nothing is dropped here. Acting on a flag means editing the registry in "
              "utils/sig_collections.py, with the reasoning written next to it.")
    else:
        print(f"\nno list flagged: none is dominated by a single gene, all behave as if they "
              f"had at least {MIN_EFFECTIVE_N_GENES} genes, and all are detected in at "
              f"least {1 - MAX_ALL_ZERO_FRACTION:.0%} of cells.")

    C.write_table(summary_tbl.round(4), "signature_concentration", coll)
    C.write_table(genes_tbl.round(6), "signature_gene_contribution", coll, index=False)

    print("\ndone.")


if __name__ == "__main__":
    main()
