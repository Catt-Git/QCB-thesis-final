#!/usr/bin/env python3
"""03_8, Route C: where the two routes agree. The main result of this stage.

A dimension-direction is called a cell state only where the two routes agree, since their
failure modes do not overlap. Verdicts: convergent (same signature family), both_different_family,
cell_only (Route A only), factor_only (Route B only; read with A_vs_null) and neither. Nothing
is promoted on a single route.

One row per live dimension and direction (dimensions from 03_6's dimension_row_order); Route
A's direction is the sign of its Spearman rho (rho > 0 on DR 7 is DR 7+).

Usage:
    export DATA_DIR=<data root>
    python convergence_tum.py                     # the scie collection, the default
    python convergence_tum.py --collection emt    # the same procedure on the EMT lists
    python convergence_tum.py --collection gavish # names the dimensions, no target axis
    python convergence_tum.py --collection gavish --all-metaprograms  # all 41, not just TNBC
"""

from __future__ import annotations

import argparse
import os
import sys

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

FDR = 0.05
N_TOP_GENES = 200        # as 03_7: the ORA table it reads

# Verdicts name which routes answered, not what that means (that depends on the rho bar).
CONVERGENT = "convergent"
BOTH_DIFFERENT = "both_different_family"
FACTOR_ONLY = "factor_only"
CELL_ONLY = "cell_only"
NEITHER = "neither"

# Labels for the figure legend and section headers.
VERDICT_LABEL = {
    CONVERGENT: "convergent",
    BOTH_DIFFERENT: "both routes, different family",
    FACTOR_ONLY: "factor only (B)",
    CELL_ONLY: "cell only (A)",
    NEITHER: "neither",
}

# p95 of Route A's null per collection (route_a_null_tum.py). Not a threshold: it only
# fills A_vs_null, i.e. whether a row below the bar is also below the noise floor.
ROUTE_A_NULL_P95 = {"scie": 0.281, "emt": 0.165, "gavish_tnbc": 0.255}

# Confounder flag thresholds on the raw-score rho; conventions, not derived.
DEPTH_FLAG = 0.30
CYCLE_FLAG = 0.30


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    return p.parse_args()


def main():
    args = parse_args()
    coll = SC.resolve(args)
    C.banner(f"03_8 - Route C, convergence: {coll.title}")
    print(f"question  {coll.question}")

    # All tables come from this collection's folder.
    rho = C.read_table("dim_signature_spearman", coll)
    # No target (gavish): no effect-size table; the convergence itself is unaffected.
    eff = C.read_table("dim_target_effect_size", coll) if coll.has_target else None
    signed = C.read_table("dim_geneset_signed_significance", coll)
    # Route B per direction from 03_7's full ORA: the signed matrix keeps only the stronger
    # side of a dimension x signature pair, so the other side would read as nothing.
    ora = pd.read_csv(C.enrichment_tsv(coll, N_TOP_GENES), sep="\t")
    ora = ora[ora["collection"] == "lab"]
    b_side = (-np.log10(ora["fdr_bh_global"].clip(lower=1e-300))).groupby(
        [ora["dimension"], ora["Term"]]).max()
    conf = C.read_table("confounders", coll)
    order = C.read_table("dimension_row_order", coll)

    dims = order.index.tolist()
    # Reported, never a filter.
    vanished = (order["vanished"].astype(bool) if "vanished" in order.columns
                else pd.Series(False, index=order.index))
    thr = -np.log10(FDR)
    # The dimension set is inherited from Route A, never recomputed.
    print(f"{len(dims)} dimensions "
          "(vanished PRUNED upstream)"
          f" x 2 directions = {2 * len(dims)} rows")
    print(f"Route A bar: |rho| >= {C.ROUTE_A_RHO_MIN};  Route B bar: global FDR < {FDR}")

    # Only the signatures both routes actually carry.
    sigs = [s for s in rho.columns if s in signed.columns]
    independent = [s for s in rho.columns if s not in signed.columns]
    if independent:
        print(f"\nRoute A also carries {', '.join(independent)}, which has no Route B "
              "counterpart by construction (it is a per-cell predictor, not a gene set).")

    rows = []
    for d in dims:
        for direction in ("+", "-"):
            s = 1.0 if direction == "+" else -1.0

            # Route A: the strongest signature association ON THIS SIDE of the axis
            a_vals = rho.loc[d, sigs].astype(float) * s
            a_best = a_vals.idxmax()
            a_rho = float(a_vals.max())

            # the same, over every Route A readout including any without a Route B counterpart
            all_vals = rho.loc[d].astype(float) * s
            a_best_any = all_vals.idxmax()
            a_rho_any = float(all_vals.max())

            # effect size of the target quadrant on this axis, oriented to the side
            if eff is None:
                smd = auroc_dir = np.nan     # no quadrant, hence no cell set to separate
            else:
                auroc = float(eff.loc[d, "auroc_target_vs_rest"])
                smd = float(eff.loc[d, "standardised_mean_difference"]) * s
                auroc_dir = auroc if direction == "+" else 1 - auroc

            # Route B: the strongest enrichment on this side
            b_vals = pd.Series({t: float(b_side.get((f"{d}{direction}", t), 0.0)) for t in sigs})
            b_best = b_vals.idxmax()
            b_neglog = float(b_vals.max())
            b_fdr = float(10 ** (-b_neglog)) if b_neglog > 0 else 1.0

            a_hit = a_rho >= C.ROUTE_A_RHO_MIN
            b_hit = b_neglog >= thr
            p95 = ROUTE_A_NULL_P95.get(coll.name)
            null_side = "unknown" if p95 is None else ("above" if a_rho >= p95 else "below")

            same_family = (coll.axis_of.get(a_best) == coll.axis_of.get(b_best)) if (a_hit and b_hit) else False
            same_signature = (a_best == b_best) if (a_hit and b_hit) else False

            if a_hit and b_hit:
                verdict = CONVERGENT if same_family else BOTH_DIFFERENT
            elif b_hit:
                verdict = FACTOR_ONLY
            elif a_hit:
                verdict = CELL_ONLY
            else:
                verdict = NEITHER

            # confounder flags, carried over from A3 for the signature being claimed
            claimed = a_best if a_hit else (b_best if b_hit else a_best)
            depth_rho = float(conf.loc[claimed, "rho_n_genes_by_counts"]) if claimed in conf.index else np.nan
            cyc = max(abs(float(conf.loc[claimed, "rho_S_score"])),
                      abs(float(conf.loc[claimed, "rho_G2M_score"]))) if claimed in conf.index else np.nan
            flags = []
            if not np.isnan(depth_rho) and abs(depth_rho) >= DEPTH_FLAG:
                flags.append("depth")
            if not np.isnan(cyc) and cyc >= CYCLE_FLAG:
                flags.append("cell_cycle")
            # Raw 03_6 rho of the claimed signature (sign = direction); a_rho is always positive here.
            claimed_rho = float(rho.loc[d, claimed]) if claimed in rho.columns else np.nan
            flags += coll.extra_flags(coll, claimed, claimed_rho)

            rows.append({
                "dimension": d, "direction": direction, "dim_direction": f"{d}{direction}",
                "A_best_signature": a_best, "A_rho": a_rho,
                "A_best_any_readout": a_best_any, "A_rho_any": a_rho_any,
                "A_auroc_target_this_side": auroc_dir,
                "A_standardised_mean_difference": smd,
                "B_best_signature": b_best, "B_neglog10_fdr": b_neglog, "B_fdr": b_fdr,
                "A_significant": a_hit, "B_significant": b_hit,
                "same_signature": same_signature, "same_family": same_family,
                "A_family": coll.axis_of.get(a_best), "B_family": coll.axis_of.get(b_best),
                "verdict": verdict,
                # Route A value against the collection's noise floor, independently of the bar.
                "A_vs_null": null_side,
                "dimension_vanished": bool(vanished[d]),
                "confounder_flags": ",".join(flags) or "none",
            })

    conv = pd.DataFrame(rows).set_index("dim_direction")
    conv = conv.loc[sorted(conv.index, key=C.dim_sort_key)]
    # Dropped rather than left as empty columns.
    auroc_col = ["A_auroc_target_this_side"]
    if not coll.has_target:
        conv = conv.drop(columns=["A_auroc_target_this_side", "A_standardised_mean_difference"])
        auroc_col = []

    C.banner("the three categories, reported separately")
    counts = conv["verdict"].value_counts()
    for v in (CONVERGENT, BOTH_DIFFERENT, FACTOR_ONLY, CELL_ONLY, NEITHER):
        print(f"  {VERDICT_LABEL[v]:32s} {int(counts.get(v, 0)):3d}")
    print("\nNo dimension is promoted on a single route. A 'convergent' row is a candidate\n"
          "cell state; the two single-route categories are candidates for the OTHER thing\n"
          "each of them can be, and are listed here for that reason, not as weaker hits.")

    C.write_table(conv, "convergence", coll)

    conv_rows = conv[conv["verdict"] == CONVERGENT].sort_values("A_rho", ascending=False)
    print(f"\nCONVERGENT ({len(conv_rows)}): both routes, same signature family")
    if len(conv_rows):
        print(conv_rows[["A_best_signature", "A_rho", *auroc_col,
                         "B_best_signature", "B_fdr", "same_signature",
                         "confounder_flags"]].to_string(float_format="%.3g"))

    b_only = conv[conv["verdict"] == FACTOR_ONLY]
    below = b_only[b_only["A_vs_null"] == "below"]
    print(f"\nFACTOR ONLY, B ({len(b_only)}): the gene program is on the axis, Route A did "
          f"not clear the bar.\n{len(below)} of them are also below the collection's own "
          f"noise floor (A_vs_null), which is\nthe subset that reads as patient-specific or "
          f"technical rather than merely weak.")
    if len(b_only):
        print(b_only[["B_best_signature", "B_fdr", "A_best_signature", "A_rho",
                      "A_vs_null", "confounder_flags"]].head(20).to_string(float_format="%.3g"))

    a_only = conv[conv["verdict"] == CELL_ONLY]
    print(f"\nCELL ONLY, A ({len(a_only)}): the cells separate, the axis does not encode the "
          "program cleanly.\nThe state may be real; this single dimension is not its description.")
    if len(a_only):
        print(a_only[["A_best_signature", "A_rho", *auroc_col,
                      "B_best_signature", "B_neglog10_fdr", "confounder_flags"]]
              .head(20).to_string(float_format="%.3g"))

    flagged = conv[(conv["verdict"] == CONVERGENT) & (conv["confounder_flags"] != "none")]
    print(f"\n{len(flagged)} of the {len(conv_rows)} convergent rows carry a confounder flag "
          "from A3 and cannot be read as clean.")

    # the project's target: both criteria of the collection on the same axis (none for gavish).
    if not coll.criteria:
        C.banner(f"{coll.name} states no target axis: it names dimensions, it does not "
                 "look for a state")
        print("The convergence table above is the result of this run. A dimension whose two\n"
              "routes independently land on the same family is a dimension with a name from\n"
              "outside this dataset; that is what a vocabulary collection is for.")
    else:
        crit_names = " AND ".join(c.label for c in coll.criteria)
        C.banner(f"the project's target for {coll.name}: {crit_names} on the same axis")

        def criterion_value(crit, d: str, side: float) -> float:
            """Strongest Route A rho meeting one criterion on one side (sign +1 high, -1 low, 0 either)."""
            names = [n for n in (crit.names or coll.by_axis(crit.axis)) if n in rho.columns]
            if not names:
                return np.nan
            v = rho.loc[d, names].astype(float) * side
            return float(v.abs().max()) if crit.sign == 0 else float((v * crit.sign).max())

        tgt = []
        for d in dims:
            for direction in ("+", "-"):
                side = 1.0 if direction == "+" else -1.0
                vals = {c.label: criterion_value(c, d, side) for c in coll.criteria}
                if all(v >= C.ROUTE_A_RHO_MIN for v in vals.values()):
                    row = {"dim_direction": f"{d}{direction}", **vals}
                    row.update({"auroc_target": conv.loc[f"{d}{direction}", "A_auroc_target_this_side"],
                                "verdict": conv.loc[f"{d}{direction}", "verdict"],
                                "flags": conv.loc[f"{d}{direction}", "confounder_flags"]})
                    tgt.append(row)
        tgt_df = pd.DataFrame(tgt)
        if len(tgt_df):
            tgt_df = tgt_df.set_index("dim_direction").sort_values("auroc_target", ascending=False)
            print(f"axes on which {crit_names} hold at the same time:")
            print(tgt_df.to_string(float_format="%.3f"))
            C.write_table(tgt_df, "target_axes", coll)
        else:
            print(f"No single axis carries both at |rho| >= {C.ROUTE_A_RHO_MIN}.\n"
                  "That is a result: the target state is an INTERSECTION of two axes in the latent\n"
                  "space rather than a direction of it, which is what the Route A target region\n"
                  "already assumed by crossing two independent scores.")

    # figures
    C.banner("figures")

    # Figure order (emt: A, B, C), as in the 03_6 and 03_7 heatmaps.
    sigs_ord = coll.order(list(sigs), for_figure=True)
    # Two panels, each with half the width cap.
    panel_w = C.fig_span(len(sigs_ord), 0.9, 3.5, cap=C.MAX_FIG_IN / 2)
    fig, axes = plt.subplots(1, 2, figsize=(2 * panel_w, C.fig_span(len(dims), 0.24, 3.5)),
                             sharey=True)
    sns.heatmap(rho.loc[dims, sigs_ord].astype(float), cmap="vlag", center=0, vmin=-0.6, vmax=0.6,
                cbar_kws={"label": "Spearman rho", "shrink": 0.4}, ax=axes[0])
    axes[0].set_title("Route A, cell-first", fontsize=10)
    vmax = float(np.nanpercentile(signed.loc[dims, sigs_ord].abs().values, 99)) or 1.0
    sns.heatmap(signed.loc[dims, sigs_ord].astype(float), cmap="vlag", center=0,
                vmin=-vmax, vmax=vmax,
                cbar_kws={"label": "signed -log10 FDR", "shrink": 0.4}, ax=axes[1])
    axes[1].set_title("Route B, factor-first", fontsize=10)
    edges = coll.block_edges(sigs_ord)
    for a in axes:
        for pos in edges:
            a.axvline(pos, color="k", lw=1.5)
        plt.setp(a.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(axes[0].get_yticklabels(), fontsize=6)
    fig.suptitle(f"Convergence, {coll.short_title}: the two routes side by side", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    C.savefig("routes_side_by_side", "03_8_convergence", coll, fig)
    plt.close(fig)

    # A strength vs B strength, one point per dimension-direction
    palette = {CONVERGENT: "#C44E52",
               BOTH_DIFFERENT: "#8172B3",
               FACTOR_ONLY: "#4C72B0",
               CELL_ONLY: "#DD8452",
               NEITHER: "0.8"}
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for v, grp in conv.groupby("verdict"):
        ax.scatter(grp["A_rho"], grp["B_neglog10_fdr"], s=34, lw=0.4, edgecolor="w",
                   c=palette.get(v, "0.5"),
                   label=f"{VERDICT_LABEL.get(v, v)} ({len(grp)})")
    ax.axvline(C.ROUTE_A_RHO_MIN, color="k", ls="--", lw=0.9)
    ax.axhline(thr, color="k", ls="--", lw=0.9)
    for lbl, r in conv[conv["verdict"] == CONVERGENT].iterrows():
        ax.annotate(lbl, (r["A_rho"], r["B_neglog10_fdr"]), fontsize=6,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(f"Route A: strongest signature association on this side (Spearman rho)")
    ax.set_ylabel(f"Route B: strongest enrichment on this side (-log10 global FDR)")
    ax.set_title(f"Convergence, {coll.short_title}: Route A vs Route B", fontsize=10)
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    sns.despine(ax=ax)
    C.savefig("convergence_scatter", "03_8_convergence", coll, fig)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
