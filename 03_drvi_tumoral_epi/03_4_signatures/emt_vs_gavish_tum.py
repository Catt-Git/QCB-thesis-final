#!/usr/bin/env python3
"""03_4: do the collaborator's EMT lists and Gavish's EMT metaprograms describe the same thing?

The emt collection (collaborator list B, a hypothesis) and the emt axis of gavish (MP12 - MP15,
a vocabulary) meet only through the cells in 03_8; this step checks whether they are simply
the same genes. The lab x Gavish Jaccard is read against two reference bands: lab x lab across
axes and Gavish EMT x Gavish EMT. A hypergeometric test complements the Jaccard. Genes are
the mapped ones, as in build_signatures_tum.py.

Usage:
    export DATA_DIR=<data root>
    python emt_vs_gavish_tum.py                        # the 4 TNBC-relevant EMT metaprograms
    python emt_vs_gavish_tum.py --all-metaprograms     # adds MP16, the glioma mesenchymal one
"""

from __future__ import annotations

import argparse
import itertools
import os
import sys

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--all-metaprograms", action="store_true",
                   help="take the EMT axis of the full 41-metaprogram registry rather than of "
                        "the TNBC-relevant subset, which adds MP16 (glioma mesenchymal)")
    return p.parse_args()


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a | b) else 0.0


def main():
    args = parse_args()
    gav = SC.GAVISH_VARIANTS[bool(args.all_metaprograms)]

    C.banner("03_4 - the collaborator's EMT lists against Gavish's EMT metaprograms")

    # inputs
    # Selected by axis, so --all-metaprograms also picks up MP16.
    mp_names = [s.name for s in gav.signatures if gav.axis_of[s.name] == "emt"]
    if not mp_names:
        sys.exit(f"[stop] the {gav.name} collection declares no signature on the 'emt' axis")

    lab_raw = C.load_signatures(SC.EMT)
    gav_raw = {n: g for n, g in C.load_signatures(gav).items() if n in mp_names}

    print(f"\ncollaborator lists : {len(lab_raw)}  ({', '.join(lab_raw)})")
    print(f"Gavish EMT metaprograms : {len(gav_raw)}  ({', '.join(gav_raw)})")

    # the universe
    # Same mapping as build_signatures_tum.py: unmapped symbols cannot count as shared.
    adata = ad.read_h5ad(C.FULL_H5AD, backed="r")
    universe = set(adata.var_names)
    n_universe = adata.n_vars
    adata.file.close()
    print(f"\nuniverse: {n_universe:,} genes of {C.FULL_H5AD.name}")

    lab = {k: set(g) & universe for k, g in lab_raw.items()}
    gavm = {k: set(g) & universe for k, g in gav_raw.items()}
    for k, g in {**lab, **gavm}.items():
        src = lab_raw.get(k, gav_raw.get(k))
        if len(g) < len(src):
            print(f"  {k:24s} {len(g):3d}/{len(src):3d} symbols mapped")

    # the matrices
    # Lab lists by axis, then the metaprograms by number.
    lab_order = SC.EMT.order(list(lab), for_figure=True)
    gav_order = [n for n in gav.order(list(gavm), for_figure=True)]
    order = lab_order + gav_order
    sets = {**lab, **gavm}

    J = pd.DataFrame(np.eye(len(order)), index=order, columns=order)
    O = pd.DataFrame(0, index=order, columns=order, dtype=int)
    for a, b in itertools.combinations(order, 2):
        J.loc[a, b] = J.loc[b, a] = jaccard(sets[a], sets[b])
        O.loc[a, b] = O.loc[b, a] = len(sets[a] & sets[b])
    for a in order:
        O.loc[a, a] = len(sets[a])

    C.write_table(J.round(4), "jaccard_emt_vs_gavish", SC.EMT)
    C.write_table(O, "shared_genes_emt_vs_gavish", SC.EMT)

    # the reference bands
    # Each band is a list of (a, b) pairs, shared by the figure and the printout.
    bands: dict[str, list[tuple[str, str]]] = {
        "lab x lab\ndifferent axis": [(a, b) for a, b in itertools.combinations(lab_order, 2)
                                      if SC.EMT.axis_of[a] != SC.EMT.axis_of[b]],
        "Gavish x Gavish\nEMT MPs only": list(itertools.combinations(gav_order, 2)),
        "lab x Gavish": [(a, b) for a in lab_order for b in gav_order],
    }

    print("\nJaccard by band (the cross-block is the question, the others are its scale)")
    summary = []
    for label, pairs in bands.items():
        v = np.array([J.loc[a, b] for a, b in pairs])
        one = label.replace("\n", ", ")
        summary.append({"band": one, "n_pairs": len(v), "median": float(np.median(v)),
                        "min": float(v.min()), "max": float(v.max())})
        print(f"  {one:34s} n={len(v):3d}   median {np.median(v):.3f}   "
              f"range {v.min():.3f} - {v.max():.3f}")
    C.write_table(pd.DataFrame(summary).set_index("band").round(4),
                  "jaccard_bands_emt_vs_gavish", SC.EMT)

    cross = np.array([J.loc[a, b] for a, b in bands["lab x Gavish"]])
    within_g = np.array([J.loc[a, b] for a, b in bands["Gavish x Gavish\nEMT MPs only"]])
    print(f"\nthe comparison that matters: lab x Gavish median {np.median(cross):.3f} against "
          f"the {np.median(within_g):.3f} that Gavish's own EMT metaprograms show each other.")
    print("Read the cross-block against THAT number, not against 1.0.")

    # the pairs, in full
    # One row per cross pair, with the shared genes spelled out.
    rows = []
    for a, b in bands["lab x Gavish"]:
        A, B = sets[a], sets[b]
        inter = sorted(A & B)
        k, n_a, n_b = len(inter), len(A), len(B)
        # P(X >= k) for k genes shared between an m-gene and an n-gene draw from n_universe.
        p = float(hypergeom.sf(k - 1, n_universe, n_a, n_b)) if k else 1.0
        rows.append({"lab_list": a, "metaprogram": b, "n_lab": n_a, "n_mp": n_b,
                     "n_shared": k, "jaccard": J.loc[a, b],
                     "overlap_coefficient": k / min(n_a, n_b),
                     "expected_by_chance": n_a * n_b / n_universe,
                     "hypergeom_p": p, "shared_genes": ";".join(inter)})
    pairs = pd.DataFrame(rows)
    # BH over the cross-block only; the reference bands are descriptive.
    o = pairs["hypergeom_p"].rank(method="first").astype(int)
    pairs["hypergeom_fdr"] = np.minimum.accumulate(
        (pairs["hypergeom_p"].sort_values(ascending=False).values
         * len(pairs) / np.arange(len(pairs), 0, -1)).clip(max=1.0)
    )[::-1][o - 1]
    pairs = pairs.sort_values("jaccard", ascending=False)
    C.write_table(pairs.round(6), "emt_vs_gavish_pairs", SC.EMT, index=False)

    print("\nthe ten strongest cross pairs:")
    print(pairs.head(10)[["lab_list", "metaprogram", "n_shared", "jaccard",
                          "expected_by_chance", "hypergeom_fdr", "shared_genes"]]
          .to_string(index=False, float_format="%.3g"))

    empty = pairs[pairs["n_shared"] == 0]
    if len(empty):
        print(f"\n{len(empty)} of the {len(pairs)} cross pairs share NO gene at all:")
        print("  " + ", ".join(f"{r.lab_list} x {r.metaprogram}" for r in empty.itertuples()))

    # figure
    # On the full matrix's scale the cross-block looks like zero, so the right panel redraws it.
    n = len(order)
    fig, (ax, axz) = plt.subplots(1, 2, figsize=(15.0, 6.6),
                                  gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.35})

    # left: the whole matrix, one scale
    off = J.where(~np.eye(n, dtype=bool))
    labelled = J.copy()
    labelled.index = labelled.columns = [f"{m} ({len(sets[m])})" for m in order]
    sns.heatmap(labelled, cmap="rocket_r", vmin=0, vmax=float(np.nanmax(off.values)),
                annot=C.annotate_cells(n, n), fmt=".2f", annot_kws={"size": 7}, square=True,
                cbar_kws={"label": "Jaccard index", "shrink": 0.6, "pad": 0.02}, ax=ax)

    # Block edges; the lab / Gavish boundary is drawn heavier.
    lab_edges = SC.EMT.block_edges(lab_order)
    for pos in lab_edges:
        ax.axhline(pos, color="black", lw=1.4)
        ax.axvline(pos, color="black", lw=1.4)
    ax.axhline(len(lab_order), color="black", lw=3)
    ax.axvline(len(lab_order), color="black", lw=3)
    ax.set_title("EMT list B vs Gavish EMT metaprograms (Jaccard)", fontsize=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)

    # right: the cross-block alone, on its own scale
    sub = J.loc[lab_order, gav_order]
    # Jaccard and shared-gene count per cell: the same Jaccard can be 6 genes or 2.
    cell = np.array([[f"{sub.loc[a, b]:.2f}\n{O.loc[a, b]:d}" for b in gav_order]
                     for a in lab_order])
    sns.heatmap(sub, cmap="rocket_r", vmin=0, vmax=float(sub.values.max()),
                annot=cell, fmt="", annot_kws={"size": 7.5}, square=False, linewidths=0.8,
                linecolor="white",
                cbar_kws={"label": "Jaccard index", "shrink": 0.75, "pad": 0.02}, ax=axz)
    axz.set_yticks(np.arange(len(lab_order)) + 0.5)
    axz.set_yticklabels([f"{a} ({len(sets[a])})" for a in lab_order], rotation=0, fontsize=8)
    axz.set_xticks(np.arange(len(gav_order)) + 0.5)
    axz.set_xticklabels([f"{b} ({len(sets[b])})" for b in gav_order], rotation=20,
                        ha="right", fontsize=8)
    for pos in SC.EMT.block_edges(lab_order):
        axz.axhline(pos, color="black", lw=1.4)
    axz.set_title("Cross-block on its own scale (Jaccard / shared genes)", fontsize=10)

    C.savefig("jaccard_emt_vs_gavish", "03_4_signatures", SC.EMT, fig)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
