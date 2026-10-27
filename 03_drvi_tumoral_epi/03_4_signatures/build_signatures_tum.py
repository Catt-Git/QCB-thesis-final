#!/usr/bin/env python3
"""03_4: ingest one collection's signature files into a .gmt, and characterise the collection.

Writes the mapped genes of each list into one .gmt (provenance in the description field)
and two tables to read before any result:

  * coverage - how much of each list is measured in this object (Route A universe, all
               genes) and in the 2,000 HVGs (Route B background). Below MIN_MAPPED_FRACTION
               the step stops: an unmapped gene is not measured, not unexpressed.
  * jaccard  - pairwise overlap: overlapping lists are not independent tests. Skipped when
               every axis has one list (emt), as there is nothing to compare.

Usage:
    export DATA_DIR=<data root>
    python build_signatures_tum.py                        # the scie collection, the default
    python build_signatures_tum.py --collection emt       # the same, on the EMT lists
    python build_signatures_tum.py --collection gavish    # the 22 TNBC-relevant metaprograms
    python build_signatures_tum.py --collection gavish --all-metaprograms   # all 41 instead
    python build_signatures_tum.py --allow-low-coverage   # report, do not stop
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
import seaborn as sns

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    p.add_argument("--allow-low-coverage", action="store_true",
                   help="report signatures below the coverage floor instead of stopping")
    return p.parse_args()


def main():
    args = parse_args()
    coll = SC.resolve(args)

    C.banner(f"03_4 - signature collection: {coll.title}")
    print(f"question    {coll.question}")
    print(f"signatures  {C.SIG_DIR}")
    print(f"object      {C.FULL_H5AD}")

    raw = C.load_signatures(coll)
    for name, genes in raw.items():
        role = "primary" if name in coll.primary_names() else "robustness"
        print(f"  {name:24s} {len(genes):5d} unique symbols  [{coll.axis_of[name]}, {role}]")

    # universes
    adata = ad.read_h5ad(C.FULL_H5AD, backed="r")
    all_genes = set(adata.var_names)
    n_all = adata.n_vars
    adata.file.close()

    hvg = ad.read_h5ad(C.HVG_H5AD, backed="r")
    hvg_genes = set(hvg.var_names)
    n_hvg = hvg.n_vars
    hvg.file.close()

    print(f"\nuniverses: {n_all:,} genes (Route A) | {n_hvg:,} HVGs (Route B ORA background)")

    # coverage
    rows = []
    mapped = {}
    for name, genes in raw.items():
        m = [g for g in genes if g in all_genes]
        mapped[name] = m
        in_hvg = [g for g in m if g in hvg_genes]
        rows.append({
            "signature": name,
            "axis": coll.axis_of[name],
            "primary": name in coll.primary_names(),
            "provenance": coll.provenance[name],
            "n_genes": len(genes),
            "n_mapped": len(m),
            "mapped_fraction": len(m) / len(genes),
            "n_in_hvg_background": len(in_hvg),
            "hvg_fraction_of_mapped": len(in_hvg) / len(m) if m else np.nan,
        })

    coverage = pd.DataFrame(rows).set_index("signature")
    coverage = coverage.loc[coll.names]

    print("\ncoverage")
    print(coverage[["axis", "primary", "n_genes", "n_mapped", "mapped_fraction",
                    "n_in_hvg_background"]].to_string(float_format="%.3f"))

    # stop checks
    low = coverage.index[coverage["mapped_fraction"] < C.MIN_MAPPED_FRACTION].tolist()
    tiny = coverage.index[coverage["n_mapped"] < C.MIN_SIGNATURE_GENES].tolist()

    if tiny:
        print(f"\n[warn] below {C.MIN_SIGNATURE_GENES} mapped genes, will be skipped "
              f"by Route A: {', '.join(tiny)}")
    if low:
        msg = (f"below the {C.MIN_MAPPED_FRACTION:.0%} mapping floor: {', '.join(low)}. "
               "Low coverage means NOT MEASURED, never not expressed - these lists carry "
               "deprecated symbols and need remapping before they mean anything.")
        if args.allow_low_coverage:
            print(f"\n[warn] {msg}")
        else:
            sys.exit(f"\n[STOP] {msg}\nRe-run with --allow-low-coverage to proceed anyway.")

    # Route B is powered by the HVG intersection, not the mapped one.
    weak_b = coverage.index[coverage["n_in_hvg_background"] < C.MIN_SIGNATURE_GENES].tolist()
    if weak_b:
        print(f"\n[warn] under {C.MIN_SIGNATURE_GENES} genes inside the 2,000-HVG ORA "
              f"background, so effectively untestable in Route B: {', '.join(weak_b)}")

    # creation of .gmt
    # Mapped genes only: the collection as actually tested.
    C.write_gmt(mapped, C.gmt_path(coll), coll.provenance)
    print(f"\n[write] {C.gmt_path(coll)}")

    C.write_table(coverage, "coverage", coll)

    order = coll.order(list(mapped), for_figure=True)
    # Jaccard asks whether lists of one axis are redundant tests; emt has one list per axis.
    if coverage["axis"].value_counts().max() < 2:
        print("\none list per axis: no Jaccard table or figure")
    else:
        jaccard_tables_and_figure(coll, mapped, order)

    # coverage barplot
    plot_coverage(coll, coverage, order, n_all)

    print("\ndone.")


def jaccard_tables_and_figure(coll, mapped: dict[str, list[str]], order: list[str]) -> None:
    """Pairwise overlap of the mapped lists: table, shared-gene counts and the heatmap."""
    names = list(mapped)
    J = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    O = pd.DataFrame(0, index=names, columns=names, dtype=int)
    for i, a in enumerate(names):
        A = set(mapped[a])
        for b in names[i + 1:]:
            B = set(mapped[b])
            inter = len(A & B)
            j = inter / len(A | B) if (A | B) else 0.0
            J.loc[a, b] = J.loc[b, a] = j
            O.loc[a, b] = O.loc[b, a] = inter
        O.loc[a, a] = len(A)

    C.write_table(J.round(4), "jaccard", coll)
    C.write_table(O, "shared_genes", coll)

    off = J.where(~np.eye(len(names), dtype=bool))
    print(f"\npairwise Jaccard: max {np.nanmax(off.values):.3f}, "
          f"median {np.nanmedian(off.values):.3f}")
    pairs = (off.stack().rename("jaccard").reset_index()
             .rename(columns={"level_0": "a", "level_1": "b"}))
    pairs = pairs[pairs["a"] < pairs["b"]].sort_values("jaccard", ascending=False)
    print("\nmost overlapping pairs (these are NOT independent readouts):")
    print(pairs.head(8).to_string(index=False, float_format="%.3f"))

    # figure
    # Mapped set size in the tick label: Jaccard alone hides very different list sizes.
    labelled = J.loc[order, order].copy()
    labelled.index = labelled.columns = [f"{n} ({len(mapped[n])})" for n in order]

    # Fixed size up to 12 lists; above, the canvas grows (capped) and the numbers come off.
    n = len(order)
    if n <= 12:
        fig_w, fig_h, annot = 8.5, 7.0, True
    else:
        fig_h = C.fig_span(n, per_item=0.55, base=3.0)
        fig_w, annot = fig_h + 1.5, C.annotate_cells(n, n)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(labelled, cmap="rocket_r", vmin=0, vmax=float(np.nanmax(off.values)),
                annot=annot, fmt=".2f", annot_kws={"size": 7}, square=True,
                cbar_kws={"label": "Jaccard index", "shrink": 0.7}, ax=ax)
    # One separator per axis boundary, with the axis name centred over its block.
    edges = coll.block_edges(order)
    blocks = list(zip([0] + edges, edges + [len(order)]))
    # Block labels go vertical on blocks under three columns; the title moves up to clear them.
    narrow = min(b - a for a, b in blocks) < 3
    ax.set_title(f"{coll.short_title}: gene overlap",
                 fontsize=10, pad=90 if narrow else 30)
    for pos in edges:
        ax.axhline(pos, color="black", lw=2)
        ax.axvline(pos, color="black", lw=2)
    for start, stop in blocks:
        ax.text((start + stop) / 2, -0.15, coll.axis_of[order[start]], ha="center",
                va="bottom", fontsize=9, weight="bold",
                rotation=90 if narrow else 0)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    C.savefig("jaccard_signature_overlap", "03_4_signatures", coll, fig)
    plt.close(fig)


def plot_coverage(coll, coverage: pd.DataFrame, order: list[str], n_all: int) -> None:
    """Mapped fraction per list, coloured by axis, with the mapping floor."""
    # 8 in up to twelve bars, wider beyond.
    fig, ax = plt.subplots(figsize=(8.0 if len(order) <= 12
                                    else C.fig_span(len(order), per_item=0.42, base=3.5), 4.5))
    cov = coverage.loc[order]
    axis_colors = dict(zip(coll.axes, C.axis_palette(len(coll.axes))))
    colors = [axis_colors[a] for a in cov["axis"]]
    ax.bar(range(len(cov)), cov["mapped_fraction"], color=colors)
    ax.axhline(C.MIN_MAPPED_FRACTION, color="crimson", ls="--", lw=1,
               label=f"floor ({C.MIN_MAPPED_FRACTION:.0%})")
    for i, (n_m, n_g) in enumerate(zip(cov["n_mapped"], cov["n_genes"])):
        ax.text(i, cov["mapped_fraction"].iloc[i] + 0.015, f"{n_m}/{n_g}",
                ha="center", fontsize=7 if len(order) <= 15 else 5, rotation=0 if len(order) <= 15 else 90)
    ax.set_xticks(range(len(cov)))
    ax.set_xticklabels(cov.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("fraction of symbols mapped")
    ax.set_ylim(0, 1.18)   # headroom for the legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=axis_colors[a]) for a in coll.axes
               if a in set(cov["axis"])]
    labels = [a for a in coll.axes if a in set(cov["axis"])]
    ax.set_title(f"{coll.short_title}: coverage ({n_all:,} genes)", fontsize=10)
    # Legend below the plot when there are many axes.
    leg_handles = handles + [ax.get_lines()[0]]
    leg_labels = labels + [f"floor ({C.MIN_MAPPED_FRACTION:.0%})"]
    if len(leg_labels) <= 6:
        ax.legend(leg_handles, leg_labels, fontsize=8)
    else:
        ax.legend(leg_handles, leg_labels, fontsize=8, loc="upper center",
                  bbox_to_anchor=(0.5, -0.42), ncol=min(8, len(leg_labels)), frameon=False)
    sns.despine(ax=ax)
    C.savefig("signature_coverage", "03_4_signatures", coll, fig)
    plt.close(fig)


if __name__ == "__main__":
    main()
