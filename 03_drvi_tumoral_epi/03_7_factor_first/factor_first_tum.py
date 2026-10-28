#!/usr/bin/env python3
"""03_7, Route B: factor-first. The dimension defines itself; prior knowledge only names it.

Each dimension-direction's top genes (OOD_combined scores of the DRVI decoder) are tested by
ORA against the collection: the discovery route, which can name dimensions nobody asked
about, but says nothing about cells. The generic naming (GO, TFs) is 03_3's factor_processes. Offline hypergeometric
tests against the 2,000-HVG training set as background, BH once across all pairs of one
collection, vanished directions pruned, rows in Route A's order.

Usage:
    export DATA_DIR=<data root>
    python factor_first_tum.py                       # the scie collection, the default
    python factor_first_tum.py --collection emt      # the same procedure on the EMT lists
    python factor_first_tum.py --collection gavish   # ORA against the TNBC metaprograms
    python factor_first_tum.py --collection gavish --all-metaprograms  # against all 41
    python factor_first_tum.py --n-top-genes 500     # a deeper list
"""

from __future__ import annotations

import argparse
import os
import sys

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gseapy as gp
import seaborn as sns
from statsmodels.stats.multitest import multipletests

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402

FDR = 0.05


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    p.add_argument("--n-top-genes", type=int, default=200,
                   help="genes per dimension-direction handed to the ORA (default 200, as in 03_3)")
    p.add_argument("--overwrite", action="store_true", help="re-run the enrichment instead of reusing the table")
    return p.parse_args()


def main():
    args = parse_args()
    coll = SC.resolve(args)
    C.banner(f"03_7 - Route B, factor-first: {coll.title}")
    print(f"question  {coll.question}")

    N_TOP = args.n_top_genes
    enrich_tsv = C.enrichment_tsv(coll, N_TOP)
    # Not collection-scoped: the top genes come off the decoder alone.
    top_tsv = C.top_genes_tsv(N_TOP)

    # ------------------------------------------------------------------ input
    embed = ad.read_h5ad(C.EMBED_H5AD)
    hvg = ad.read_h5ad(C.HVG_H5AD, backed="r")
    gene_names = hvg.var_names.copy()
    hvg.file.close()

    assert embed.varm[f"{C.SCORE_KEY}_positive"].shape[1] == len(gene_names), \
        "the embedding and the DRVI input disagree on the gene axis"

    n_van = C.n_vanished(embed)
    print(f"{embed.n_obs:,} cells x {embed.n_vars} latent dimensions "
          f"({n_van} flagged vanished in var['vanished'], PRUNED per direction below)")
    print(f"{len(gene_names):,} HVGs, the DRVI training feature set and the ORA background")

    # -------------------------------------------- per-direction top gene lists
    scores_df = C.interpretability_scores(embed, gene_names, key=C.SCORE_KEY)
    n_dir = scores_df.shape[1]
    print(f"\n{scores_df.shape[0]:,} genes x {n_dir} dimension-directions "
          f"({2 * embed.n_vars - n_dir} of 2 x {embed.n_vars} dropped as vanished IN THAT "
          f"DIRECTION; the BH denominator below follows)")
    print("Both directions of each dimension are tested separately: DRVI can encode two\n"
          "distinct concepts on the two sides of one axis, and pooling them cancels the\n"
          "two programs against each other.")

    top_genes = {d: scores_df[d].sort_values(ascending=False).head(N_TOP).index.tolist()
                 for d in scores_df.columns}
    pd.DataFrame(top_genes, index=pd.RangeIndex(1, N_TOP + 1, name="rank")).to_csv(top_tsv, sep="\t")
    print(f"[write] {top_tsv}")

    # ------------------------------------------------------------- gene sets
    custom = C.read_gmt(C.gmt_path(coll))
    # ORA can only ever test the part of a signature that lives in the background.
    custom_bg = {k: [g for g in v if g in set(gene_names)] for k, v in custom.items()}
    print("\nsignature sizes inside the ORA background (the tested form):")
    for k, v in custom_bg.items():
        flag = "  <- under the floor, effectively untestable" if len(v) < C.MIN_SIGNATURE_GENES else ""
        print(f"  {k:24s} {len(custom[k]):5d} mapped -> {len(v):4d} in background{flag}")

    # One library, the requested collection: BH runs over its pairs alone.
    libraries = {"lab": custom_bg}

    background = list(gene_names)

    # ------------------------------------------------------------- enrichment
    if enrich_tsv.exists() and not args.overwrite:
        print(f"\n[read] {enrich_tsv} (--overwrite to re-run)")
        long = pd.read_csv(enrich_tsv, sep="\t")
    else:
        C.banner("ORA, offline hypergeometric against the 2,000-HVG background")
        records = []
        for i, (dim, genes) in enumerate(top_genes.items(), start=1):
            for lib, sets in libraries.items():
                try:
                    res = gp.enrich(gene_list=genes, gene_sets=sets,
                                    background=background, outdir=None)
                except Exception as exc:
                    print(f"  [{i:>3}/{len(top_genes)}] {dim:>8} {lib:9s} FAILED: {exc}")
                    continue
                # gseapy returns an empty list, not a DataFrame, when nothing overlaps.
                res_df = getattr(res, "results", None)
                if res_df is None or len(res_df) == 0:
                    continue
                df = pd.DataFrame(res_df).copy()
                df.insert(0, "collection", lib)
                df.insert(0, "dimension", dim)
                records.append(df)
            if i % 20 == 0 or i == len(top_genes):
                print(f"  [{i:>3}/{len(top_genes)}] queried")

        long = pd.concat(records, ignore_index=True)
        long.to_csv(enrich_tsv, sep="\t", index=False)
        print(f"[write] {enrich_tsv}")

    long["direction"] = long["dimension"].str[-1]
    long["dim"] = long["dimension"].str[:-1].str.strip()

    # ------------------------------------------------------------------ BH
    # Once across all pairs tested; pairs with no overlap are absent from gseapy's output but
    # belong in the denominator, so they are added back as p = 1.
    C.banner("Benjamini-Hochberg across all dimension-direction / gene-set pairs")
    n_terms = sum(len(s) for s in libraries.values())
    n_pairs = len(top_genes) * n_terms
    print(f"{len(top_genes)} directions x {n_terms} gene sets = {n_pairs:,} pairs tested; "
          f"{len(long):,} returned a non-empty overlap")

    pvals = np.concatenate([long["P-value"].values, np.ones(n_pairs - len(long))])
    rej, padj, _, _ = multipletests(pvals, alpha=FDR, method="fdr_bh")
    long["fdr_bh_global"] = padj[:len(long)]
    long["significant"] = long["fdr_bh_global"] < FDR
    print(f"{int(long['significant'].sum()):,} pairs significant at global FDR < {FDR} "
          f"(gseapy's own per-query adjustment would have called "
          f"{int((long['Adjusted P-value'] < FDR).sum()):,})")

    print(f"\nNOTE for the Methods: the {len(custom_bg)} signatures of this collection are NOT\n"
          "independent - see the Jaccard matrix of 03_4, whose block structure is exactly\n"
          "this. BH assumes independence or positive dependence; the correction here is\n"
          f"therefore conservative in count but must NOT be read as {len(custom_bg)} independent\n"
          f"tests of {len(custom_bg)} independent hypotheses.")

    long.to_csv(enrich_tsv, sep="\t", index=False)

    # ------------------------------------------------- signed matrix, Route A order
    order_tbl = C.read_table("dimension_row_order", coll)
    row_order = order_tbl.index.tolist()
    print(f"\nrow order taken from Route A: {len(row_order)} dimensions")

    sig_cols = coll.order(list(custom_bg))
    # The table keeps registry order, the figure the collection's figure order (emt: A, B, C).
    fig_cols = coll.order(list(custom_bg), for_figure=True)
    signed = pd.DataFrame(0.0, index=row_order, columns=sig_cols)
    sign_of = pd.DataFrame("", index=row_order, columns=sig_cols)

    h = long[long["collection"] == "lab"].copy()
    h["neglog"] = -np.log10(h["fdr_bh_global"].clip(lower=1e-300))
    for (dim, term), grp in h.groupby(["dim", "Term"], observed=True):
        if dim not in signed.index or term not in signed.columns:
            continue
        best = grp.loc[grp["neglog"].idxmax()]
        signed.loc[dim, term] = best["neglog"] * (1 if best["direction"] == "+" else -1)
        sign_of.loc[dim, term] = best["direction"]

    C.write_table(signed.round(4), "dim_geneset_signed_significance", coll)
    C.write_table(long[long["significant"]].drop(columns=["Genes"], errors="ignore"),
                  "factor_first_significant", coll, index=False)

    thr = -np.log10(FDR)
    n_hit = (signed.abs() >= thr).sum().sort_values(ascending=False)
    print("\nsignificant dimension-directions per signature (global FDR)")
    print(n_hit.to_string())

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(figsize=(C.fig_span(len(fig_cols), 1.0, 4.0),
                                    C.fig_span(len(row_order), 0.24, 3.0)))
    vmax = float(np.nanpercentile(signed.abs().values, 99)) or 1.0
    sns.heatmap(signed[fig_cols], cmap="vlag", center=0, vmin=-vmax, vmax=vmax,
                cbar_kws={"label": "signed -log10 FDR  (+ = positive direction)", "shrink": 0.4},
                ax=ax)
    for pos in coll.block_edges(fig_cols):
        ax.axvline(pos, color="k", lw=1.5)
    ax.set_title(f"Factor-first, {coll.short_title}: DRVI dimensions x gene sets", fontsize=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    C.savefig("dim_geneset_signed_heatmap", "03_7_factor_first", coll, fig)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
