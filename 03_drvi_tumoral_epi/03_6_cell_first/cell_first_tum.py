#!/usr/bin/env python3
"""03_6, Route A: cell-first. Prior knowledge defines the state, DRVI is the coordinate system.

Each collection's readouts are scored per cell (all genes, unintegrated object), standardised
within cohort, and correlated with every live DRVI dimension. The only route that assigns
cells to states. The target region comes from the collection (utils/sig_collections.py):
scie takes stem-high / immunogenic-low, emt the hybrid list high (E x M co-expression shown beside it);
gavish defines none, so A5 and its figures are skipped.

Named checks on the target set: the cell cycle (recomputed inside G1, for every collection),
depth for an immune-low group (scie) and doublet / ambient RNA for a mesenchymal-high group
(emt). The region is defined once per plane and the stability of the cell set across the
definitions is itself a result. A1 - A5 never read the embedding; only A6 does.

Usage:
    export DATA_DIR=<data root>
    python cell_first_tum.py                              # the scie collection, the default
    python cell_first_tum.py --collection emt             # the same procedure on the EMT lists
    python cell_first_tum.py --collection gavish          # scores only: no target region
    python cell_first_tum.py --collection gavish --all-metaprograms  # all 41, not the TNBC 22
    python cell_first_tum.py --high-q 0.80 --low-q 0.20   # a stricter target region
    python cell_first_tum.py --overwrite                  # re-score instead of reusing the csv
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
import scanpy as sc
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402

N_BINS = 25                 # sc.tl.score_genes default

# Cohort only: inside the malignant compartment the groupings are states, and standardising
# within them would remove the contrast measured. They are reported as covariates instead.
GROUPBY = ["cohort"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    C.add_cutoff_arguments(p)
    p.add_argument("--overwrite", action="store_true", help="re-score even if the per-cell csv exists")
    return p.parse_args()



# Scoring
def score_signatures(adata, sets: dict[str, list[str]]) -> tuple[pd.DataFrame, list[str]]:
    """One sc.tl.score_genes call per signature, ctrl_size matched to its length, seeded."""
    scores, skipped = {}, []
    for name, genes in sets.items():
        if len(genes) < C.MIN_SIGNATURE_GENES:
            print(f"  {name:24s} SKIPPED, only {len(genes)} mapped genes "
                  f"(floor {C.MIN_SIGNATURE_GENES})")
            skipped.append(name)
            continue
        sc.tl.score_genes(
            adata, gene_list=genes, ctrl_size=len(genes), n_bins=N_BINS,
            random_state=C.SEED, use_raw=False, score_name=f"score_{name}",
        )
        v = adata.obs[f"score_{name}"].astype(float)
        scores[name] = v.values
        print(f"  {name:24s} {len(genes):5d} genes  mean {v.mean():+.4f}  sd {v.std():.4f}")
    return pd.DataFrame(scores, index=adata.obs_names), skipped


def standardise_within(df: pd.DataFrame, strata: pd.DataFrame) -> pd.DataFrame:
    """Z-score each column within the strata of GROUPBY; a stratum with no spread gets 0."""
    g = df.groupby([strata[k].values for k in GROUPBY], observed=True)
    z = g.transform(lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else s * 0.0)
    return z.fillna(0.0).add_prefix("z_")



# Quadrants
# The region rules live in signature_common, shared with later steps.

def axis_label(coll, name: str) -> str:
    """The axis a readout sits on, marking the ones that do not come from a gene list (CytoTRACE2)."""
    if name in coll.extra_readouts:
        return f"{coll.extra_readouts[name]} (independent)"
    return coll.axis_of.get(name, "")


def jaccard(a: pd.Series, b: pd.Series) -> float:
    inter = int((a & b).sum())
    union = int((a | b).sum())
    return inter / union if union else np.nan


# Share of the definitions a cell needs for the one-panel plane (4 of 6 for scie); stricter
# than the 03_6 consensus, which takes ceil(n/2).
MEAN_PLANE_VOTE_FRACTION = 2 / 3


def plane_mean_figure(z: pd.DataFrame, qdf: pd.DataFrame, planes, cut, coll, idx) -> None:
    """One panel: the x readouts of the voting planes averaged, red = called by most definitions.

    Only when every voting plane shares its y readout and rules (scie: stemness high x the
    immune readout low); otherwise the mean of the x axes is not one axis.
    """
    if len(planes) < 3 or len({(pl.y, pl.y_rule, pl.x_rule) for pl in planes}) != 1 \
            or planes[0].y is None:
        return
    y_name, y_rule, x_rule = planes[0].y, planes[0].y_rule, planes[0].x_rule
    n = len(planes)
    k = int(np.ceil(MEAN_PLANE_VOTE_FRACTION * n))
    called = qdf[[pl.label for pl in planes]].sum(axis=1) >= k
    # Mean of the within-stratum z-scores, not re-standardised: 0 is still the stratum mean.
    zx = z[[f"z_{pl.x}" for pl in planes]].mean(axis=1)
    zy = z[f"z_{y_name}"]
    print(f"\none-panel plane: {int(called.sum()):,} cells called by >= {k} of {n} definitions")

    x, y, tgt = zx.values[idx], zy.values[idx], called.values[idx]
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.scatter(x[~tgt], y[~tgt], s=1.5, c="0.78", lw=0, rasterized=True)
    ax.scatter(x[tgt], y[tgt], s=1.5, c="#C44E52", lw=0, rasterized=True)
    q = {"high": cut.high_q, "low": cut.low_q}
    ax.axhline(zy.quantile(q[y_rule]), color="k", ls="--", lw=0.8)
    # The red set is not cut on the mean: this line is only a reference.
    ax.axvline(zx.quantile(q[x_rule]), color="0.45", ls=":", lw=0.8)
    ax.set_xlabel(f"mean z of the {n} {coll.axis_of.get(planes[0].x, 'x')} readouts "
                  f"({x_rule}) (within {' x '.join(GROUPBY)})", fontsize=8)
    ax.set_ylabel(f"z {y_name} ({y_rule})", fontsize=8)
    ax.set_title(f"{int(called.sum()):,} cells called by at least {k} of {n} definitions",
                 fontsize=9)
    ax.tick_params(labelsize=7)
    sns.despine(ax=ax)
    fig.text(0.01, 0.005, f"dashed: q{q[y_rule]:.2f} of {y_name}; dotted: q{q[x_rule]:.2f} of "
             f"the mean, a reference (red is the vote, not this cut). Readouts: "
             f"{', '.join(pl.x for pl in planes)}.", fontsize=6, color="0.35", wrap=True)
    fig.suptitle(f"Cell-first, {coll.short_title}: target region, one panel", fontsize=10)
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    C.savefig(f"{coll.plane_figure}_mean", "03_6_cell_first", coll, fig)
    plt.close(fig)


def main():
    args = parse_args()
    coll = SC.resolve(args)
    C.banner(f"03_6 - Route A, cell-first: {coll.title}")
    print(f"question  {coll.question}")
    print(f"space     DRVI ({C.RUN_ID})")

    def savefig_close(name, step, coll_, fig):
        path = C.savefig(name, step, coll_, fig)
        plt.close(fig)
        return path

    sets = C.read_gmt(C.gmt_path(coll))
    print(f"{len(sets)} signatures from {C.gmt_path(coll)}")

    # input
    adata = ad.read_h5ad(C.FULL_H5AD)
    C.attach_leiden(adata.obs)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes (all genes, unintegrated)")
    # Must be the malignant object of 03_2.
    assert adata.obs[C.LABEL_KEY].astype(str).nunique() >= 1, "no cell_type column"
    assert (adata.obs["cnv_status"].astype(str) != "not_tested").all(), \
        "not_tested cells in the object: this is not the 03_2 subset"

    # `.X` must be log-normalised: score_genes on raw counts is a different statistic.
    xmax = float(adata.X.max())
    assert xmax < 50, (f".X looks like counts (max {xmax:.1f}), not log-normalised data. "
                       "Route A must score on the scran log1p matrix.")
    print(f".X max {xmax:.3f}, log1p in uns: {'log1p' in adata.uns}  -> log-normalised")

    for k in ("S_score", "G2M_score", "phase", "n_genes_by_counts", "pct_counts_mt"):
        assert k in adata.obs, f"missing obs key: {k}"
    print("cell cycle scores reused from 03_2, not recomputed and not regressed out")

    # scores
    scores_csv = C.scores_csv(coll)
    if scores_csv.exists() and not args.overwrite:
        print(f"\n[read] {scores_csv} (--overwrite to re-score)")
        allsc = pd.read_csv(scores_csv, index_col=0)
        assert (allsc.index == adata.obs_names).all(), "score csv is not in the object's cell order"
        raw = allsc[[c for c in allsc.columns if not c.startswith("z_")]]
        raw.columns = [c.replace("score_", "") for c in raw.columns]
        # Derived readouts are rebuilt from the z-scores below.
        raw = raw.drop(columns=[d.name for d in coll.derived], errors="ignore")
        skipped = [n for n in sets if n not in raw.columns and n not in coll.extra_readouts]
    else:
        print("\nscoring (sc.tl.score_genes, use_raw=False, ctrl_size=len(signature), "
              f"n_bins={N_BINS}, random_state={C.SEED})")
        raw, skipped = score_signatures(adata, sets)

    # CytoTRACE2
    # Only scie declares it (a stemness readout).
    ct2_col = None
    if "CytoTRACE2" in coll.extra_readouts and C.CYTOTRACE_CSV.exists():
        ct2 = pd.read_csv(C.CYTOTRACE_CSV, index_col=0).reindex(adata.obs_names)
        if "CytoTRACE2_Score" in ct2:
            raw["CytoTRACE2"] = ct2["CytoTRACE2_Score"].values
            ct2_col = "CytoTRACE2"
            print(f"\n[read] {C.CYTOTRACE_CSV}: CytoTRACE2_Score joined as an "
                  "extra stemness readout, independent of the signature lists")
    elif "CytoTRACE2" in coll.extra_readouts:
        print(f"\n[warn] {C.CYTOTRACE_CSV} not found - CytoTRACE2 is left out of the "
              "quadrant definitions.\n"
              "       Run 03_5_cytotrace2/cytotrace2_tum.py in the cytotrace2-py env.")

    # standardisation
    # Before the confounders: the derived readouts are contrasts of z-scores.
    C.banner(f"A3 - standardisation within ({', '.join(GROUPBY)})")
    z = standardise_within(raw, adata.obs)
    n_strata = adata.obs.groupby(GROUPBY, observed=True).size()
    print(f"{len(n_strata)} strata, smallest {n_strata.min()} cells, largest {n_strata.max():,}")
    print(f"singleton strata (z forced to 0): {int((n_strata == 1).sum())}")

    # derived readouts
    # z(plus) - z(minus) on the standardised scores, so both halves share a scale.
    for d in coll.derived:
        missing = [n for n in (d.plus, d.minus) if f"z_{n}" not in z.columns]
        if missing:
            print(f"  {d.name:24s} SKIPPED, missing {', '.join(missing)}")
            continue
        v = z[f"z_{d.plus}"] - z[f"z_{d.minus}"]
        z[f"z_{d.name}"] = v
        raw[d.name] = v
        print(f"  {d.name:24s} {d.description}")

    readout_axes = {n: axis_label(coll, n) for n in raw.columns}
    for ax_name in coll.axes:
        on_axis = [n for n in raw.columns if coll.axis_of.get(n) == ax_name]
        if on_axis:
            print(f"\n{ax_name:12s} readouts: {', '.join(on_axis)}")

    # A vocabulary collection (gavish) has no target: everything on a quadrant is skipped,
    # not emptied, so no table reports a target of zero cells.
    if coll.has_target:
        planes = coll.planes(coll, list(raw.columns))
        assert any(pl.in_consensus for pl in planes), (
            f"the {coll.name} collection defines no plane on the readouts that were "
            "scored - nothing can be called a target region")
        for pl in planes:
            for k in pl.readouts:
                assert k in raw.columns, f"plane '{pl.label}' needs readout {k}, which was not scored"
        print(f"\n{len(planes)} definitions of the target region ({coll.target_label}):")
        for pl in planes:
            print(f"  {pl.label:24s} {pl.describe()}")
    else:
        planes = []
        print(f"\n{coll.name} defines no target region: it is a vocabulary, not a hypothesis.")
        print("SKIPPED: A5 (the region, its consensus, its stability, its per-patient and\n"
              "         per-group breakdowns), the named risks computed on it, the per-dimension\n"
              "         target effect size, and the plane and stability figures.")
        print("PRODUCED: the per-cell scores, the confounder table, and the dimensions x\n"
              "          readouts correlations - everything that needs no region to mean\n"
              "          something. Route B and Route C run unchanged.")

    # confounder table
    C.banner("A4 - confounders")
    conf_keys = ["n_genes_by_counts", "pct_counts_mt", "S_score", "G2M_score"]
    rows = []
    for name in raw.columns:
        r = {"readout": name, "axis": readout_axes[name]}
        for k in conf_keys:
            rho, p = spearmanr(raw[name].values, adata.obs[k].values)
            r[f"rho_{k}"] = rho
            r[f"p_{k}"] = p
        rows.append(r)
    conf = pd.DataFrame(rows).set_index("readout")
    print(conf[[f"rho_{k}" for k in conf_keys]].to_string(float_format="%+.3f"))
    C.write_table(conf, "confounders", coll)

    cc_max = conf[["rho_S_score", "rho_G2M_score"]].abs().max(axis=1)
    print(f"\nstrongest cell-cycle coupling: {cc_max.idxmax()} (|rho| = {cc_max.max():.3f})")
    print(f"strongest depth coupling     : {conf['rho_n_genes_by_counts'].abs().idxmax()} "
          f"(|rho| = {conf['rho_n_genes_by_counts'].abs().max():.3f})")

    out = pd.concat([raw.add_prefix("score_"), z], axis=1)
    out.to_csv(scores_csv)
    print(f"\n[write] {scores_csv}")

    if coll.has_target:
        # A5 - target region
        C.banner(f"A5 - the target region: {coll.target_label}\n"
                 f"(high >= q{args.high_q:.2f}, low <= q{args.low_q:.2f}, "
                 f"mid = q{args.mid_lo_q:.2f} - q{args.mid_hi_q:.2f})")

        labels = [pl.label for pl in planes]
        # Only these vote on the consensus; the others are compared with it.
        voting = [pl.label for pl in planes if pl.in_consensus]
        cut = C.Cutoffs.from_args(args)
        quads = {pl.label: C.define_target(z, pl, cut) for pl in planes}

        qdf = pd.DataFrame(quads)
        sizes = qdf.sum().rename("n_cells").to_frame()
        sizes["pct_of_compartment"] = 100 * sizes["n_cells"] / adata.n_obs
        print(f"\ntarget region ({coll.target_label}), one definition per plane")
        print(sizes.to_string(float_format="%.2f"))

        # Stability across definitions; for emt, the hybrid list against the co-expression corner.
        stab = None
        if len(labels) > 1:
            stab = pd.DataFrame(index=labels, columns=labels, dtype=float)
            for a in labels:
                for b in labels:
                    stab.loc[a, b] = jaccard(qdf[a], qdf[b])
            print("\nstability of the cell set across definitions (Jaccard)")
            print(stab.to_string(float_format="%.3f"))
            off = stab.where(~np.eye(len(labels), dtype=bool))
            print(f"median pairwise Jaccard: {np.nanmedian(off.values):.3f}  "
                  f"(range {np.nanmin(off.values):.3f} - {np.nanmax(off.values):.3f})")
            C.write_table(stab.round(4), "quadrant_stability", coll)
        else:
            print("\none definition only: no stability table")

        n_defs, consensus = C.consensus_vote(qdf[voting])
        print(f"\ncalled by >=1 definition: {int((n_defs >= 1).sum()):,} cells; "
              f"by a majority: {int(consensus.sum()):,}; by all {len(voting)}: "
              f"{int((n_defs == len(voting)).sum()):,}")
        for lab in [lb for lb in labels if lb not in voting]:
            both = int((qdf[lab] & consensus).sum())
            print(f"{lab} (no vote): {int(qdf[lab].sum()):,} cells, {both:,} of them in the target")

        # A state present in one patient is a patient effect until shown otherwise.
        per_pat = qdf.copy()
        per_pat["cohort"] = adata.obs["cohort"].values
        pp = per_pat.groupby("cohort", observed=True).sum()
        pp["n_cells_in_cohort"] = adata.obs["cohort"].value_counts().reindex(pp.index).values
        pp["consensus"] = pd.Series(consensus.values, index=adata.obs["cohort"].values).groupby(level=0).sum()
        pp["pct_consensus"] = 100 * pp["consensus"] / pp["n_cells_in_cohort"]
        print("\nper-patient size of the consensus quadrant "
              "(a state in one patient is a patient effect until shown otherwise)")
        print(pp[["n_cells_in_cohort", "consensus", "pct_consensus"]]
              .sort_values("pct_consensus", ascending=False).to_string(float_format="%.2f"))
        print(f"\npatients with at least one consensus cell: "
              f"{int((pp['consensus'] > 0).sum())} / {len(pp)}")
        C.write_table(pp, "quadrant_per_patient", coll)

        # Leiden is a covariate of the target set: reported, never standardised away.
        for key in [k for k in C.grouping_keys(adata.obs) if k == C.LEIDEN_KEY]:
            per_g = pd.DataFrame({key: adata.obs[key].astype(str).values,
                                  "consensus": consensus.values}).groupby(key, observed=True).agg(
                n_cells=("consensus", "size"), n_target=("consensus", "sum"))
            per_g["pct"] = 100 * per_g["n_target"] / per_g["n_cells"]
            print(f"\nconsensus quadrant by {key}")
            print(per_g.sort_values("pct", ascending=False).to_string(float_format="%.2f"))
            C.write_table(per_g, f"quadrant_per_{key}", coll)

        # the named risks (A4 cont.)
        C.banner(f"A4 - the named risks of this collection: {', '.join(coll.risks)}")

        # 1. is the target set just the cycling one? Runs for every collection.
        phase = adata.obs["phase"].astype(str).values
        comp = pd.crosstab(pd.Series(phase, name="phase"), consensus.values,
                           normalize="columns") * 100
        comp.columns = ["rest", "target"]
        print("phase composition, target quadrant vs the rest (%)")
        print(comp.to_string(float_format="%.2f"))

        g1 = phase == "G1"
        print(f"\nrecomputing the target region inside G1 alone ({g1.sum():,} cells), "
              "i.e. with the cycle held out")
        # Cutoffs recomputed within G1: what the definition calls without the cycling cells.
        z_g1 = z[g1]
        g1_quads = {pl.label: C.define_target(z_g1, pl, cut) for pl in planes if pl.in_consensus}
        g1df = pd.DataFrame(g1_quads)
        _, g1_consensus = C.consensus_vote(g1df)
        overlap = jaccard(consensus[g1], g1_consensus)
        print(f"consensus target restricted to G1: {int(consensus[g1].sum()):,} cells")
        print(f"consensus target recomputed within G1: {int(g1_consensus.sum()):,} cells")
        print(f"Jaccard between the two: {overlap:.3f}")
        print("A high Jaccard means the state is not an artefact of the cycle; a low one means the\n"
              "target was largely 'cycling' and the readout does not survive the check.")

        cc_rows = [{"check": "phase_pct_G1_target", "value": float(comp.loc["G1", "target"]) if "G1" in comp.index else np.nan},
                   {"check": "phase_pct_G1_rest", "value": float(comp.loc["G1", "rest"]) if "G1" in comp.index else np.nan},
                   {"check": "n_target_all_phases", "value": float(consensus.sum())},
                   {"check": "n_target_within_G1_recomputed", "value": float(g1_consensus.sum())},
                   {"check": "jaccard_target_vs_G1_recomputed", "value": float(overlap)}]

        # 2. is the low end of the primary axis just shallow sequencing? (`risks` contains "depth")
        if "depth" in coll.risks and coll.depth_risk_readout:
            r = coll.depth_risk_readout
            low_grp = z[f"z_{r}"] <= z[f"z_{r}"].quantile(args.low_q)
            depth = adata.obs["n_genes_by_counts"].values
            u, pu = mannwhitneyu(depth[low_grp.values], depth[~low_grp.values], alternative="two-sided")
            med_lo, med_hi = np.median(depth[low_grp.values]), np.median(depth[~low_grp.values])
            auc_depth = roc_auc_score(low_grp.values, -depth)
            print(f"\n{r}-low ({int(low_grp.sum()):,} cells) vs the rest, n_genes_by_counts:")
            print(f"  median {med_lo:,.0f} vs {med_hi:,.0f}  (Mann-Whitney p = {pu:.3g})")
            print(f"  AUROC of 'shallower' predicting {r}-low: {auc_depth:.3f}")
            print("  0.5 means depth does not explain the group; well above it means the evasive\n"
                  "  group is largely the low-complexity group and the finding is technical.")
            cc_rows += [{"check": "median_depth_immunogenic_low", "value": float(med_lo)},
                        {"check": "median_depth_rest", "value": float(med_hi)},
                        {"check": "mannwhitney_p_depth", "value": float(pu)},
                        {"check": "auroc_depth_predicts_immunogenic_low", "value": float(auc_depth)}]

        # 3. is the high end of the mesenchymal axis just ambient RNA or a doublet? Aneuploid cells
        # cannot be fibroblasts, but ambient fibroblast RNA survives subsetting (doublet_score: 01_2).
        if "ambient" in coll.risks and coll.ambient_risk_axis:
            on_axis = [n for n in raw.columns if coll.axis_of.get(n) == coll.ambient_risk_axis]
            dbl = adata.obs["doublet_score"].astype(float).values
            print(f"\n{coll.ambient_risk_axis}-high vs the rest, doublet_score (Scrublet, via 03_2):")
            for r in on_axis:
                high = (z[f"z_{r}"] >= z[f"z_{r}"].quantile(args.high_q)).values
                auc = roc_auc_score(high, dbl)
                rho_d, p_d = spearmanr(raw[r].values, dbl)
                med_hi_d, med_lo_d = np.median(dbl[high]), np.median(dbl[~high])
                print(f"  {r:24s} AUROC {auc:.3f}   median {med_hi_d:.4f} vs {med_lo_d:.4f}   "
                      f"rho(score, doublet_score) {rho_d:+.3f}")
                cc_rows += [{"check": f"auroc_doublet_predicts_{r}_high", "value": float(auc)},
                            {"check": f"rho_doublet_score_{r}", "value": float(rho_d)}]
            if "predicted_doublet" in adata.obs:
                pred = adata.obs["predicted_doublet"].astype(bool).values
                print(f"  predicted doublets in the consensus target: "
                      f"{int(pred[consensus.values].sum()):,} / {int(consensus.sum()):,} "
                      f"({100 * pred[consensus.values].mean():.2f}%) vs "
                      f"{100 * pred[~consensus.values].mean():.2f}% in the rest")
                cc_rows += [{"check": "pct_predicted_doublet_target",
                             "value": float(100 * pred[consensus.values].mean())},
                            {"check": "pct_predicted_doublet_rest",
                             "value": float(100 * pred[~consensus.values].mean())}]
            print("  AUROC near 0.5 means the mesenchymal signal is not the soup; well above it\n"
                  "  means the EMT readout is measuring contamination and nothing here is a state.")

        C.write_table(pd.DataFrame(cc_rows).set_index("check"), "confounder_checks", coll, index=True)

    # A6 - dimensions x signatures
    C.banner("A6 - the DRVI dimensions vs the standardised scores")
    if not C.EMBED_H5AD.exists():
        raise SystemExit(
            f"{C.EMBED_H5AD} not found.\nRun 03_3_drvi_run/run_drvi_tum.py (or drvi_tum_64.ipynb) first.")
    embed = ad.read_h5ad(C.EMBED_H5AD)
    assert (embed.obs_names == adata.obs_names).all(), \
        "the embedding and the all-genes object do not hold the same cells in the same order"
    dims = C.analysis_dimensions(embed)
    n_van = C.n_vanished(embed)
    print(f"{embed.n_vars} dimensions, {len(dims)} used")
    # 03_7 reads the row order off the table below.
    dropped = [d for d in embed.var["title"] if d not in set(dims)]
    print(f"{n_van} flagged vanished in var['vanished'] and "
          f"{len(dropped)} DROPPED -> {', '.join(dropped)}")

    L = pd.DataFrame(np.asarray(embed.X), index=embed.obs_names, columns=embed.var["title"].values)[dims]

    readouts = list(raw.columns)
    rho = pd.DataFrame(index=dims, columns=readouts, dtype=float)
    for r in readouts:
        rr, _ = spearmanr(L.values, z[f"z_{r}"].values)
        rho[r] = rr[:-1, -1] if rr.ndim else np.nan
    # A dimension with no spread gives NaN: set to 0, no association.
    rho = rho.fillna(0.0)
    C.write_table(rho.round(4), "dim_signature_spearman", coll)
    print(f"\nstrongest |rho| per readout:")
    for r in readouts:
        print(f"  {r:24s} {rho[r].abs().idxmax():>8s}  rho = {rho.loc[rho[r].abs().idxmax(), r]:+.3f}")

    # Effect size of the consensus quadrant on each dimension; only with a target.
    if coll.has_target:
        auroc = pd.Series({d: roc_auc_score(consensus.values, L[d].values) for d in dims},
                          name="auroc_target_vs_rest")
        # sd == 0 only on a dead dimension: 0, not inf.
        smd = pd.Series({d: ((L[d][consensus.values].mean() - L[d][~consensus.values].mean())
                             / L[d].std(ddof=0) if L[d].std(ddof=0) > 0 else 0.0)
                         for d in dims}, name="standardised_mean_difference")
        eff = pd.concat([auroc, smd], axis=1)
        eff["abs_auroc_from_half"] = (eff["auroc_target_vs_rest"] - 0.5).abs()
        eff = eff.sort_values("abs_auroc_from_half", ascending=False)
        print("\neffect size of the consensus target quadrant on each dimension (top 10)")
        print(eff.head(10).to_string(float_format="%+.3f"))
        C.write_table(eff, "dim_target_effect_size", coll)

    # The row order every later heatmap uses, Route B included.
    order = pd.Series(dims, name="dimension").to_frame().assign(
        drvi_order=embed.var.set_index("title").loc[dims, "order"].values,
        vanished=embed.var.set_index("title").loc[dims, "vanished"].astype(bool).values)
    C.write_table(order.set_index("dimension"), "dimension_row_order", coll)

    # figures
    C.banner("figures")

    # The two figures of the target region.
    if coll.has_target:
        # the plane, one panel per definition of the target region
        def cut_lines(a, v, rule, vertical: bool):
            """Draw the cutoff(s) of one rule; 'mid' has two, so it reads as a band."""
            draw = a.axvline if vertical else a.axhline
            qs = {"high": [args.high_q], "low": [args.low_q],
                  "mid": [args.mid_lo_q, args.mid_hi_q]}[rule]
            for q in qs:
                draw(v.quantile(q), color="k", ls="--", lw=0.8)

        ncol = min(4, len(planes))
        nrow = int(np.ceil(len(planes) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.8 * nrow), squeeze=False)
        rng = np.random.default_rng(C.SEED)
        idx = rng.choice(adata.n_obs, size=min(20000, adata.n_obs), replace=False)
        for a, pl in zip(axes.ravel(), planes):
            zx = z[f"z_{pl.x}"]
            called = qdf[pl.label]
            if pl.y:
                zy = z[f"z_{pl.y}"]
                x, y = zx.values[idx], zy.values[idx]
                tgt = called.values[idx]
                a.scatter(x[~tgt], y[~tgt], s=1.5, c="0.78", lw=0, rasterized=True)
                a.scatter(x[tgt], y[tgt], s=1.5, c="#C44E52", lw=0, rasterized=True)
                cut_lines(a, zy, pl.y_rule, vertical=False)
                a.set_ylabel(f"z {pl.y} ({pl.y_rule})", fontsize=8)
            else:
                # One readout: every cell, the called ones stacked in red.
                bins = np.histogram_bin_edges(zx.values, bins=80)
                a.hist([zx.values[~called.values], zx.values[called.values]], bins=bins,
                       stacked=True, color=["0.78", "#C44E52"], lw=0)
                a.set_ylabel("cells", fontsize=8)
            cut_lines(a, zx, pl.x_rule, vertical=True)
            note = "" if pl.in_consensus else f", {int((called & consensus).sum()):,} in the target"
            a.set_title(f"{pl.label}\n{int(called.sum()):,} cells{note}", fontsize=9)
            a.set_xlabel(f"z {pl.x} ({pl.x_rule}) (within {' x '.join(GROUPBY)})", fontsize=8)
            a.tick_params(labelsize=7)
            sns.despine(ax=a)
        for a in axes.ravel()[len(planes):]:
            a.axis("off")
        fig.suptitle(f"Cell-first, {coll.short_title}: target region (red)", fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        savefig_close(coll.plane_figure, "03_6_cell_first", coll, fig)
        # The same plane in one panel, beside the per-definition one.
        plane_mean_figure(z, qdf, [pl for pl in planes if pl.in_consensus], cut, coll, idx)

        # quadrant stability; two definitions are one number, printed and in the table
        if stab is not None and len(labels) > 2:
            fig, ax = plt.subplots(figsize=(6, 5))
            sns.heatmap(stab.astype(float), cmap="rocket_r", vmin=0, vmax=1, annot=True, fmt=".2f",
                        annot_kws={"size": 7}, square=True,
                        cbar_kws={"label": "Jaccard of the called cell set", "shrink": 0.7}, ax=ax)
            ax.set_title(f"Cell-first, {coll.short_title}: target stability", fontsize=10)
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
            plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
            savefig_close("quadrant_stability", "03_6_cell_first", coll, fig)

    # dimensions x signatures: the sign of rho gives the direction (rho > 0 on DR 7 is DR 7+).
    # Derived readouts (EMT_SCORE_B) are contrasts of the columns drawn: kept for 03_8, not drawn.
    derived_names = {d.name for d in coll.derived}
    col_order = coll.order([c for c in rho.columns if c not in derived_names], for_figure=True)
    fig, ax = plt.subplots(figsize=(C.fig_span(len(col_order), 1.0, 4.0),
                                    C.fig_span(len(dims), 0.24, 3.0)))
    sns.heatmap(rho[col_order].astype(float), cmap="vlag", center=0, vmin=-0.6, vmax=0.6,
                cbar_kws={"label": "Spearman rho (dimension vs within-stratum z-score)\n"
                                   "sign = direction: rho > 0 is DR n+, rho < 0 is DR n-",
                          "shrink": 0.4}, ax=ax)
    ax.set_title(f"Cell-first, {coll.short_title}: DRVI dimensions x signatures", fontsize=10)
    for pos in coll.block_edges(col_order):
        ax.axvline(pos, color="k", lw=1.5)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    C.savefig("dim_signature_heatmap", "03_6_cell_first", coll, fig)
    plt.close(fig)

    if skipped:
        print(f"\n[warn] signatures skipped for having under {C.MIN_SIGNATURE_GENES} "
              f"mapped genes: {', '.join(skipped)}")
    print("\ndone.")


if __name__ == "__main__":
    main()
