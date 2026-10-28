#!/usr/bin/env python3
"""03_8, the DRVI-Leiden clusters: what is in them, and which directions they carry.

Reads the 03_3 embedding and its Leiden partitions on the live DRs, leiden_smi_tum.py's SMI
matrix, STATE_PICKS and the naming tables. The clusters are cut from the space they are
described in, so everything here is descriptive and circular. One mode per run:

  --states-de      the STATE_PICKS clusters: DE vs rest on all genes, per-patient consistency,
                   composition, markers vs decoder top genes, ORA
  --states-dr      every cluster: which live directions are enriched (AUROC), what each
                   direction is, one-row summary per cluster (CLUSTER_READING). Reads the DE

The DE and DR functions take any partition: 03_8b runs them on the max(DR) groups.

Outputs: tables/03_8_convergence/state_*_drvi_leiden_<res>_<run>, and
figures/03_8_convergence/state_de_drvi_leiden_<res>/, state_dr_enrichment_drvi_leiden_<res>/.

Usage:
    export DATA_DIR=<data root>
    python3 leiden_states_tum.py --states-de
    python3 leiden_states_tum.py --states-dr
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import warnings

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "utils"))
sys.path.insert(0, HERE)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402
from dr_naming_tum import STATE_PICKS  # noqa: E402
from leiden_smi_tum import load_partitions  # noqa: E402

CS = C.CS

STEP = "03_8_convergence"

FDR = 0.05                       # as 03_7, and applied the same way

N_TOP_GENES = 200                # as 03_7: the list depth is shared by both arms

STATE_DE_MIN_LFC = 0.5      # log2 fold change a characterisation marker has to clear

STATE_DE_MIN_PCT = 0.10     # and the share of the cluster's cells that must express it

STATE_DE_MIN_CELLS = 20     # patient x (cluster | rest) cells for the per-patient check

STATE_DE_N_CHECK = 50       # top markers per cluster checked patient by patient

STATE_DR_AUROC = 0.80       # a direction counts as ENRICHED in a cluster from here

STATE_DR_TOP = 3            # directions per cluster carried into the summary

# Curated reading of each cluster at Leiden 1.5, from the --states-de / --states-dr outputs of
# 2026-10-27 (draft, to be checked by eye): the one factor that best explains it. Patient shares
# as in FACTOR_STYLE; object median 1,408 genes per cell. Re-read if the partition changes.
#   factor: program | single_patient | few_patients | depth | contamination | low_quality | none
CLUSTER_READING = {
    "0": ("single_patient", "interferon / MHC (HLA-B, LY6E, HLA-DRA)"),
    "1": ("depth", "depth axis DR 1+, fibroblast ambient (SPARC, DCN, COL1A2), 653 genes"),
    "2": ("none", "ribosomal-rich, low genes (RPL18A, RPS12), 566 genes; DR 3+"),
    "3": ("few_patients", "CD24, PCSK1N, NUPR1; DR 21+ senescence"),
    "4": ("few_patients", "secretory (WFDC2, SLPI, DSC2)"),
    "5": ("program", "cell cycle G1/S (STMN1, H2AFZ, TUBB), DR 5+"),
    "6": ("program", "hypoxia (BNIP3, P4HA1, FAM162A), DR 12-"),
    "7": ("single_patient", "MHC-I / II (HLA-A/B/C, CD74, HLA-DRA)"),
    "8": ("low_quality", "nuclear lncRNA (MALAT1, NEAT1, KCNQ1OT1), 745 genes; DR 6-"),
    "9": ("program", "luminal secretory, MP23 SECRETED 2 (MUCL1, AZGP1, APOD), DR 17-"),
    "10": ("single_patient", "VIM, S100B, IGFBP5; depth axis DR 4+"),
    "11": ("single_patient", "baseline sample; VIM, BMP7, S100B, 4,408 genes"),
    "12": ("low_quality", "immediate-early / dissociation stress (ZFP36, FOS, IER2)"),
    "13": ("none", "housekeeping / ribosomal (FTH1, EEF1A1), 719 genes"),
    "14": ("single_patient", "interferon / S100 (CXCL10, WARS, SERPING1, S100A14)"),
    "15": ("single_patient", "RT+PD1 sample; CRIP2, SMIM22, LYPD3, WFDC2"),
    "16": ("contamination", "ambient haemoglobin (HBB, HBA1, HBA2), 419 genes; DR 33+"),
    "17": ("low_quality", "nuclear lncRNA + stress (NEAT1, MALAT1, PPP1R15A, FOSB)"),
    "18": ("contamination", "myeloid (C1QA-C, APOE, TYROBP) on depth axis DR 4+"),
    "19": ("few_patients", "basal / epithelial (KRT5, DSP, SPINT2), 3,677 genes"),
    "20": ("none", "ribosomal + luminal (RPS4X, XBP1, TFF3, AZGP1); DR 30-"),
    "21": ("few_patients", "LIM_STEM / EMT_1 axis (DR 20-, DR 11+); CSRP2, IER3"),
    "22": ("low_quality", "mitochondrial transcripts (MT-CYB, MT-CO3, MT-ATP6)"),
    "23": ("few_patients", "secretory (FDCSP, AZGP1, WFDC2), DR 35+"),
    "24": ("program", "cell cycle G2/M (CDC20, CCNB1, UBE2C), DR 2+"),
    "25": ("single_patient", "secretory (PRR27, CSN3, ODAM)"),
    "26": ("single_patient", "secretory (PRR27, CSN3, ODAM, STATH); DR 11-"),
    "27": ("few_patients", "S100A8/A9 senescence (S100A8, S100A9, S100A7), DR 32+"),
    "28": ("program", "MHC-II / IFN-gamma (CD74, HLA-DRA, HLA-B), DR 15-"),
    "29": ("contamination", "immune (PTPRC, CCL5, SRGN, CXCR4), doublets / ambient; DR 19-"),
    "30": ("single_patient", "baseline sample; CRABP2, RBP1, PDLIM1, 3,971 genes; DR 41-"),
    "31": ("single_patient", "stromal-like (MGP, TAGLN, SPARC); DR 20-"),
    "32": ("single_patient", "POMC, NFATC4, AR, HOXB2, 4,638 genes; DR 44+"),
    "33": ("single_patient", "secretory (CSN3, PRR27, TSPAN13)"),
}

# Slots 1-5 of the categorical palette; greys for the absence of an explanation.
FACTOR_STYLE = {
    "program": ("#2a78d6", "biological program, across patients"),
    "single_patient": ("#eb6834", "single patient (>= 80% of cells)"),
    "few_patients": ("#f5b48f", "few patients (top patient 50-80% of cells)"),
    "depth": ("#1baf7a", "sequencing depth"),
    "contamination": ("#eda100", "contamination (ambient / immune)"),
    "low_quality": ("#e87ba4", "low quality / dissociation stress"),
    "none": ("#9a9a9a", "no clear program"),
    "too_small": ("#cfcfcf", "too small to read (< 50 cells)"),
}

INK = "#222222"


def states_paths(step: str = STEP):
    return C.fig_dir(step), C.table_dir(step)


def live_directions(embed) -> pd.DataFrame:
    """Cells x live DR directions, each clipped at zero: `DR k+` = max(x, 0), `DR k-` =
    max(-x, 0). Vanished directions dropped, as in the SMI and the AUROC."""
    X = np.asarray(embed.X)
    titles = embed.var["title"].astype(str).to_numpy()
    cols, vals = [], []
    for j, t in enumerate(titles):
        if not bool(embed.var["vanished_positive_direction"].iloc[j]):
            cols.append(f"{t}+"); vals.append(X[:, j].clip(min=0))
        if not bool(embed.var["vanished_negative_direction"].iloc[j]):
            cols.append(f"{t}-"); vals.append((-X[:, j]).clip(min=0))
    return pd.DataFrame(np.column_stack(vals), index=embed.obs_names.astype(str), columns=cols)


def _direction_order(d: str):
    """`DR 3+` < `DR 3-` < `DR 10+`; anything that is not a direction goes last."""
    try:
        return (int(d.split()[1][:-1]), d[-1] == "-")
    except (IndexError, ValueError):
        return (10 ** 6, False)


def leiden_partition(args, drvi_run: str) -> dict:
    """The groups --states-de and --states-dr read: the Leiden clusters at --smi-resolution.

    Returns key (output name), lab (cell -> group), groups, picks (characterised by DE),
    tested (DE run on them), links (group -> its directions) and the curated reading.
    """
    parts = load_partitions(drvi_run)
    res = args.smi_resolution
    if str(res) not in parts.columns:
        sys.exit(f"resolution {res} is not one of 03_3's: have {list(parts.columns)}")
    lab = parts[str(res)].astype(int).astype(str)
    groups = [str(c) for c in sorted(lab.astype(int).unique())]
    return {"key": f"drvi_leiden_{res}", "lab": lab, "unit": "cluster", "groups": groups,
            "picks": list(STATE_PICKS), "tested": groups, "links": STATE_PICKS,
            "desc": f"Leiden on the live DRs, resolution {res}",
            "reading": CLUSTER_READING, "what": "DRVI-Leiden cluster"}


def drvi_states_de(hvg, args, drvi_run: str, P: dict | None = None, step: str = STEP):
    """Characterise the picked clusters: what is in them, and is it their DR's program.

    1. DE vs the rest, Wilcoxon, all genes expressed in >= 1% of cells;
    2. per-patient consistency of the top markers (a cluster can be one patient);
    3. composition: patient, treatment, response, phase, depth / mt / doublet / CNV;
    4. markers (HVG, top N) against each direction's top decoder genes, hypergeometric;
    5. ORA against scie, emt and gavish_tnbc, BH per library.
    `P` is the partition (default the Leiden clusters); outputs go to `step`.
    """
    from scipy.stats import hypergeom

    P = P or leiden_partition(args, drvi_run)
    key, picks, links, unit = P["key"], P["picks"], P["links"], P["unit"]
    fdir, tdir = states_paths(step)
    ddir = fdir / f"state_de_{key}"
    ddir.mkdir(exist_ok=True)
    C.banner(f"DRVI states: characterising {len(picks)} {unit}s of {key}")

    # the cells
    full = ad.read_h5ad(C.FULL_H5AD)
    full.obs_names = full.obs_names.astype(str)
    embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad", backed="r")
    meta = embed.obs.copy()
    meta.index = meta.index.astype(str)
    embed.file.close()
    meta = meta.reindex(full.obs_names)
    lab = P["lab"].reindex(full.obs_names).astype(str)
    full.obs[key] = pd.Categorical(lab.values, categories=P["groups"])
    sc.pp.filter_genes(full, min_cells=int(0.01 * full.n_obs))
    print(f"    {full.n_obs:,} cells x {full.n_vars:,} genes (expressed in >= 1% of cells)",
          flush=True)

    # 1. DE
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(full, groupby=key, groups=P["tested"], method="wilcoxon",
                                use_raw=False, reference="rest", pts=True, key_added="_de")
    de = sc.get.rank_genes_groups_df(full, group=None, key="_de")
    de["group"] = de["group"].astype(str)
    de = de.rename(columns={"group": "cluster", "names": "gene"})
    print(f"    DE: {len(P['tested'])} {unit}s vs rest in {time.time() - t0:.0f}s", flush=True)
    sig = de[(de["pvals_adj"] < 0.05) & (de["logfoldchanges"].abs() > STATE_DE_MIN_LFC)]
    up = sig[(sig["logfoldchanges"] > 0) & (sig["pct_nz_group"] >= STATE_DE_MIN_PCT)]
    up = up.sort_values(["cluster", "scores"], ascending=[True, False])
    sig.to_csv(tdir / f"state_de_{key}_{drvi_run}.csv.gz", index=False)
    print(f"[table] {tdir / f'state_de_{key}_{drvi_run}.csv.gz'}  (all clusters, "
          f"|lfc| > {STATE_DE_MIN_LFC}, FDR < 0.05)")

    # 2. per patient
    cohort = meta["cohort"].astype(str)
    rows = []
    for cl in picks:
        genes = up[up["cluster"] == cl]["gene"].head(STATE_DE_N_CHECK).tolist()
        if not genes:
            continue
        X = full[:, genes].X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        inside = (lab == cl).values
        n_ok = 0
        deltas = []
        for pat in cohort.unique():
            pm = (cohort == pat).values
            a, b = pm & inside, pm & ~inside
            if a.sum() < STATE_DE_MIN_CELLS or b.sum() < STATE_DE_MIN_CELLS:
                continue
            n_ok += 1
            deltas.append(X[a].mean(0) - X[b].mean(0))
        deltas = np.array(deltas) if deltas else np.zeros((0, len(genes)))
        for j, g in enumerate(genes):
            rows.append({"cluster": cl, "gene": g, "rank": j + 1, "n_patients_tested": n_ok,
                         "n_patients_up": int((deltas[:, j] > 0).sum()) if n_ok else 0})
    cons = pd.DataFrame(rows)
    cons["share_patients_up"] = cons["n_patients_up"] / cons["n_patients_tested"].replace(0, np.nan)
    marker_tbl = up[up["cluster"].isin(picks)].merge(cons, on=["cluster", "gene"], how="left")
    marker_tbl.to_csv(tdir / f"state_markers_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'state_markers_{key}_{drvi_run}.csv'}")

    # 3. composition
    comp = []
    for cl in picks:
        m = (lab == cl).values
        pat = cohort[m].value_counts(normalize=True)
        row = {"cluster": cl, "directions": ", ".join(links[cl]), "n_cells": int(m.sum()),
               "n_patients": int((cohort[m].value_counts() >= 1).sum()),
               "top_patient": pat.index[0], "top_patient_share": float(pat.iloc[0]),
               "patients_for_80pct": int((pat.cumsum() < 0.8).sum() + 1)}
        for col in ("treatment", "response", "phase"):
            for k, v in meta.loc[m, col].astype(str).value_counts(normalize=True).items():
                row[f"{col}_{k}"] = float(v)
        for col in ("n_genes_by_counts", "pct_counts_mt", "doublet_score", "cnv_score"):
            row[f"median_{col}"] = float(meta.loc[m, col].median())
            row[f"median_{col}_all"] = float(meta[col].median())
        comp.append(row)
    comp = pd.DataFrame(comp).set_index("cluster")
    comp.to_csv(tdir / f"state_composition_{key}_{drvi_run}.csv")
    print(f"[table] {tdir / f'state_composition_{key}_{drvi_run}.csv'}")

    # 4. link to the DRs
    hvg_genes = list(hvg.var_names)
    hv = set(hvg_genes)
    n_top = args.n_top_genes
    dec = pd.read_csv(C.top_genes_tsv(n_top), sep="\t", index_col=0)
    emb_var = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad", backed="r").var
    live = set(emb_var.loc[~emb_var["vanished_positive_direction"].astype(bool), "title"] + "+") | \
        set(emb_var.loc[~emb_var["vanished_negative_direction"].astype(bool), "title"] + "-")
    dec = dec[[c for c in dec.columns if c in live]]
    lists = {cl: [g for g in up[up["cluster"] == cl]["gene"] if g in hv][:n_top] for cl in picks}
    ov = []
    for cl, genes in lists.items():
        gs = set(genes)
        for d in dec.columns:
            dg = set(dec[d].dropna())
            k = len(gs & dg)
            p = float(hypergeom.sf(k - 1, len(hv), len(dg), len(gs))) if k else 1.0
            ov.append({"cluster": cl, "direction": d, "n_markers_hvg": len(gs), "overlap": k,
                       "p": p, "linked": d in links[cl]})
    ov = pd.DataFrame(ov)
    ov["fdr_bh"] = multipletests(ov["p"], method="fdr_bh")[1]
    ov["rank_in_cluster"] = ov.groupby("cluster")["p"].rank(method="min")
    ov.to_csv(tdir / f"state_marker_vs_decoder_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'state_marker_vs_decoder_{key}_{drvi_run}.csv'}")

    # 5. ORA
    libs = {c.name: C.read_gmt(C.gmt_path(c))
            for c in (SC.get("scie"), SC.get("emt"), SC.resolve(argparse.Namespace(
                collection="gavish", all_metaprograms=False)))}
    ora = []
    for lib, sets in libs.items():
        sets_bg = {k: [g for g in v if g in hv] for k, v in sets.items()}
        sets_bg = {k: v for k, v in sets_bg.items() if len(v) >= C.MIN_SIGNATURE_GENES}
        long = run_ora({cl: g for cl, g in lists.items() if g}, sets_bg, hvg_genes, lib, key)
        long, _ = apply_bh(long, sum(1 for g in lists.values() if g), len(sets_bg))
        long["library"] = lib
        ora.append(long)
    ora = pd.concat(ora, ignore_index=True)
    ora = ora.drop(columns=["Genes"], errors="ignore")
    ora.to_csv(tdir / f"state_ora_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'state_ora_{key}_{drvi_run}.csv'}")

    # report
    C.banner("per cluster")
    for cl in picks:
        c = comp.loc[cl]
        mk = marker_tbl[marker_tbl["cluster"] == cl]
        o = ov[ov["cluster"] == cl].sort_values("p")
        link = o[o["linked"]]
        sig_ora = ora[(ora["gene_list"] == cl) & ora["significant"]].sort_values("fdr_bh")
        print(f"\n{unit} {cl}  ({c['n_cells']:,} cells, {c['directions']})")
        print(f"  patients: {c['n_patients']}, top {c['top_patient']} "
              f"{c['top_patient_share']:.0%}, {c['patients_for_80pct']} make 80%")
        tr = ", ".join(f"{a} {c.get(f'treatment_{a}', 0):.0%}" for a in ("BASE", "PD1", "RTPD1"))
        print(f"  treatment: {tr}")
        print(f"  median genes {c['median_n_genes_by_counts']:.0f} (all "
              f"{c['median_n_genes_by_counts_all']:.0f}), mt {c['median_pct_counts_mt']:.1f}% "
              f"(all {c['median_pct_counts_mt_all']:.1f}%), doublet "
              f"{c['median_doublet_score']:.2f} (all {c['median_doublet_score_all']:.2f})")
        print(f"  {len(mk)} up markers; top 15: " + ", ".join(
            f"{r.gene}({r.logfoldchanges:.1f}"
            + (f";{int(r.n_patients_up)}/{int(r.n_patients_tested)}p" if r.n_patients_tested == r.n_patients_tested else "")
            + ")" for r in mk.head(15).itertuples()))
        print("  decoder overlap, linked: " + "; ".join(
            f"{r.direction} {r.overlap}/{r.n_markers_hvg} p={r.p:.1e} (rank {int(r.rank_in_cluster)})"
            for r in link.itertuples()))
        print("  decoder overlap, best 3: " + "; ".join(
            f"{r.direction} {r.overlap} p={r.p:.1e}" for r in o.head(3).itertuples()))
        print("  ORA (FDR<0.05): " + ("; ".join(
            f"{r.library}:{r.Term} {r.fdr_bh:.1e}" for r in sig_ora.head(6).itertuples())
            or "none"))

    # figures
    # (a) dot plot, top 6 markers per picked cluster (2 for the max(DR) groups), the rest pooled as 'other'.
    n_dot = 6 if len(picks) <= 12 else 2
    top6 = {cl: [g for g in marker_tbl[marker_tbl["cluster"] == cl]["gene"].head(n_dot)]
            for cl in picks}
    grp = lab.where(lab.isin(picks), "other")
    full.obs["_pick"] = pd.Categorical(grp.values, categories=picks + ["other"])
    dp = sc.pl.dotplot(full, var_names={f"{k}" if unit == "group" else f"cl {k}": v
                                        for k, v in top6.items() if v},
                       groupby="_pick", standard_scale="var", show=False, return_fig=True,
                       title=f"top markers (one-vs-rest, all genes) of the picked {unit}s, {key}")
    dp.savefig(ddir / f"dotplot_markers_{key}_{drvi_run}.png", dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"[fig] {ddir / f'dotplot_markers_{key}_{drvi_run}.png'}")

    # (b) composition: patient and treatment shares per picked cluster, all cells for scale
    fig, axes = plt.subplots(1, 2, figsize=(13, (0.45 if len(picks) <= 12 else 0.26)
                                            * (len(picks) + 1) + 1.8),
                             gridspec_kw={"width_ratios": [3, 1.2]})
    rows_ = picks + ["all"]
    pats = cohort.value_counts().index.tolist()
    pal = dict(zip(pats, plt.get_cmap("tab20")(np.linspace(0, 1, len(pats)))))
    for ax, col, order, colours in (
            (axes[0], cohort, pats, pal),
            (axes[1], meta["treatment"].astype(str), ["BASE", "PD1", "RTPD1"],
             {"BASE": "#2a78d6", "PD1": "#eb6834", "RTPD1": "#1baf7a"})):
        for i, cl in enumerate(rows_):
            m = np.ones(len(lab), bool) if cl == "all" else (lab == cl).values
            share = col[m].value_counts(normalize=True).reindex(order).fillna(0)
            left = 0.0
            for k, v in share.items():
                ax.barh(i, v, left=left, color=colours[k], height=0.7,
                        edgecolor="white", lw=0.5, label=k if i == 0 else None)
                left += v
        ax.set_yticks(range(len(rows_)), [f"{unit} {c}" if c != "all" else "all cells"
                                          for c in rows_], fontsize=8 if len(rows_) <= 13 else 6)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("share of cells", fontsize=8)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes[0].set_title("patient", fontsize=9)
    axes[1].set_title("treatment", fontsize=9)
    axes[0].legend(fontsize=6, ncol=1, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    axes[1].legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    fig.suptitle(f"Convergence, DRVI states: composition of the picked {unit}s", fontsize=10)
    fig.tight_layout()
    fig.savefig(ddir / f"composition_{key}_{drvi_run}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {ddir / f'composition_{key}_{drvi_run}.png'}")


def _auroc(values: np.ndarray, inside: np.ndarray) -> float:
    """AUROC of `values` separating `inside` from the rest, via ranks (ties averaged)."""
    from scipy.stats import rankdata
    n1, n0 = int(inside.sum()), int((~inside).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(values)
    return float((r[inside].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def direction_interpretation(drvi_run: str, live: list[str]) -> pd.DataFrame:
    """One row per live direction: what the phase already says about it, nothing recomputed.

    Name (03_8), depth axis, Route B hits, 03_3's GO term, top decoder genes;
    `interpretation` is the first that exists.
    """
    run = drvi_run
    naming = pd.read_csv(C.table_path("dr_naming", "consensus", run), comment="#", index_col=0)
    progs = pd.read_csv(C.table_path("program_dimensions", "consensus", run), comment="#")
    named = {}
    for r in progs.itertuples():
        for d in str(r.dimensions).split(","):
            named.setdefault(d.strip(), []).append(r.program)
    route_b = {}
    for coll in ("scie", "emt", "gavish_tnbc"):
        cv = pd.read_csv(C.table_path("convergence", coll, run), comment="#", index_col=0)
        for d, r in cv[cv["B_significant"]].iterrows():
            route_b.setdefault(d, []).append((float(r["B_fdr"]),
                                              f"{r['B_best_signature']} (FDR {r['B_fdr']:.0e}, "
                                              f"A rho {r['A_rho']:.2f})"))
    # 03_3's generic naming: the Enrichr GO term, else g:Profiler's.
    fp = pd.read_csv(C.factor_processes_csv(), index_col="factor")
    go_best = fp["enrichr_term"].fillna(fp["gprofiler_term"]).dropna()
    dec = pd.read_csv(C.top_genes_tsv(N_TOP_GENES), sep="\t", index_col=0)

    rows = []
    for d in live:
        dim = d[:-1].strip()
        cov = naming.loc[dim] if dim in naming.index else None
        depth = cov is not None and abs(float(cov["rho_log_umi"])) >= 0.30
        nm = named.get(d, [])
        rb = [t for _, t in sorted(route_b.get(d, []))]     # strongest Route B hit first
        go = go_best.get(d, "")
        genes = ", ".join(dec[d].head(10)) if d in dec.columns else ""
        if nm:
            interp = "named: " + ", ".join(nm)
        elif depth:
            interp = f"depth axis (dimension rho log UMI {float(cov['rho_log_umi']):+.2f})"
        elif rb:
            interp = "Route B only: " + rb[0].split(" (")[0]
        elif go:
            interp = "GO only: " + go
        else:
            interp = "unnamed - read the genes"
        rows.append({
            "direction": d, "interpretation": interp,
            "named_programs": "; ".join(nm), "route_b_hits": "; ".join(rb),
            "go_term_03_3": go, "top10_decoder_genes": genes,
            "dim_rho_log_umi": float(cov["rho_log_umi"]) if cov is not None else np.nan,
            "dim_rho_S": float(cov["rho_S_score"]) if cov is not None else np.nan,
            "dim_rho_G2M": float(cov["rho_G2M_score"]) if cov is not None else np.nan,
            "dim_rho_pct_mt": float(cov["rho_pct_mt"]) if cov is not None else np.nan,
        })
    return pd.DataFrame(rows).set_index("direction")


def drvi_states_dr(args, drvi_run: str, P: dict | None = None, step: str = STEP):
    """DR enrichment per cluster, for every cluster of one resolution.

    AUROC of each direction's signed coordinate, cluster vs rest (the two sides sum to 1); for the
    top STATE_DR_TOP directions also the median within-patient AUROC. Circular: what it adds is
    which program, confounder or patient each cluster is made of.
    """
    P = P or leiden_partition(args, drvi_run)
    key, unit = P["key"], P["unit"]
    fdir, tdir = states_paths(step)
    odir = fdir / f"state_dr_enrichment_{key}"
    odir.mkdir(exist_ok=True)
    de_path = tdir / f"state_de_{key}_{drvi_run}.csv.gz"
    if not de_path.exists():
        sys.exit(f"missing {de_path}: run --states-de first")
    de = pd.read_csv(de_path)
    de["cluster"] = de["cluster"].astype(str)
    C.banner(f"DRVI states: DR enrichment and interpretation, every {unit} of {key}")

    embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad")
    embed.obs_names = embed.obs_names.astype(str)
    meta = embed.obs
    lab = P["lab"].reindex(embed.obs_names).astype(str)
    clusters = list(P["tested"])
    X = np.asarray(embed.X)
    titles = embed.var["title"].astype(str).to_numpy()
    cols, vals = [], []
    for j, t in enumerate(titles):
        if not bool(embed.var["vanished_positive_direction"].iloc[j]):
            cols.append(f"{t}+"); vals.append(X[:, j])
        if not bool(embed.var["vanished_negative_direction"].iloc[j]):
            cols.append(f"{t}-"); vals.append(-X[:, j])
    V = np.column_stack(vals)
    print(f"    {len(clusters)} {unit}s x {len(cols)} live directions", flush=True)

    # AUROC, pooled
    from scipy.stats import rankdata
    ranks = np.column_stack([rankdata(V[:, k]) for k in range(V.shape[1])])
    au = pd.DataFrame(index=clusters, columns=cols, dtype=float)
    smd = pd.DataFrame(index=clusters, columns=cols, dtype=float)
    sd = V.std(0)
    for cl in clusters:
        m = (lab == cl).values
        n1, n0 = m.sum(), (~m).sum()
        au.loc[cl] = (ranks[m].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)
        smd.loc[cl] = (V[m].mean(0) - V[~m].mean(0)) / np.where(sd > 0, sd, 1)
    smi_path = C.table_dir("03_8_convergence") / f"smi_matrix_{key}_{drvi_run}.csv"
    smi = pd.read_csv(smi_path, index_col=0) if smi_path.exists() else None

    interp = direction_interpretation(drvi_run, cols)
    interp.to_csv(tdir / f"state_dr_interpretation_{key}_{drvi_run}.csv")
    print(f"[table] {tdir / f'state_dr_interpretation_{key}_{drvi_run}.csv'}")

    # long table
    cohort = meta["cohort"].astype(str).values
    long = []
    for cl in clusters:
        m = (lab == cl).values
        order = au.loc[cl].sort_values(ascending=False)
        for rank, (d, a) in enumerate(order.items(), start=1):
            row = {"cluster": cl, "direction": d, "rank": rank, "auroc": float(a),
                   "smd": float(smd.loc[cl, d]),
                   "smi": float(smi.loc[d, cl]) if smi is not None and d in smi.index
                   and cl in smi.columns else np.nan,
                   "enriched": bool(a >= STATE_DR_AUROC)}
            if rank <= STATE_DR_TOP:
                k = cols.index(d)
                within = []
                for pat in np.unique(cohort):
                    pm = cohort == pat
                    if (pm & m).sum() >= STATE_DE_MIN_CELLS and (pm & ~m).sum() >= STATE_DE_MIN_CELLS:
                        within.append(_auroc(V[pm, k], m[pm]))
                row["n_patients_within"] = int(len(within))
                row["median_auroc_within"] = float(np.median(within)) if within else np.nan
            long.append(row)
    long = pd.DataFrame(long).merge(interp[["interpretation"]], left_on="direction",
                                    right_index=True, how="left")
    long.to_csv(tdir / f"state_dr_enrichment_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'state_dr_enrichment_{key}_{drvi_run}.csv'}")

    # one row per cluster
    up = de[(de["logfoldchanges"] > STATE_DE_MIN_LFC) & (de["pct_nz_group"] >= STATE_DE_MIN_PCT)]
    up = up.sort_values(["cluster", "scores"], ascending=[True, False])
    summ = []
    for cl in clusters:
        m = (lab == cl).values
        pat = pd.Series(cohort[m]).value_counts(normalize=True)
        tr = meta.loc[m, "treatment"].astype(str).value_counts(normalize=True)
        top = long[(long["cluster"] == cl) & (long["rank"] <= STATE_DR_TOP)]
        summ.append({
            "cluster": cl, "n_cells": int(m.sum()),
            "n_patients_ge20": int((pd.Series(cohort[m]).value_counts() >= 20).sum()),
            "top_patient": pat.index[0], "top_patient_share": float(pat.iloc[0]),
            "BASE": float(tr.get("BASE", 0)), "PD1": float(tr.get("PD1", 0)),
            "RTPD1": float(tr.get("RTPD1", 0)),
            "median_genes": float(meta.loc[m, "n_genes_by_counts"].median()),
            "median_pct_mt": float(meta.loc[m, "pct_counts_mt"].median()),
            "median_doublet": float(meta.loc[m, "doublet_score"].median()),
            "n_enriched_directions": int(((long["cluster"] == cl) & long["enriched"]).sum()),
            "top_directions": "; ".join(
                f"{r.direction} (AUROC {r.auroc:.2f}, within {r.median_auroc_within:.2f} "
                f"in {int(r.n_patients_within)}p)" for r in top.itertuples()),
            "top_direction_reading": "; ".join(f"{r.direction}: {r.interpretation}"
                                               for r in top.itertuples()),
            "top_markers": ", ".join(up[up["cluster"] == cl]["gene"].head(10)),
        })
    summ = pd.DataFrame(summ).set_index("cluster")
    summ.to_csv(tdir / f"state_cluster_summary_{key}_{drvi_run}.csv")
    print(f"[table] {tdir / f'state_cluster_summary_{key}_{drvi_run}.csv'}")

    C.banner("per cluster")
    med_g = float(meta["n_genes_by_counts"].median())
    for cl, r in summ.iterrows():
        print(f"\n{unit} {cl}: {r.n_cells:,} cells; {r.top_patient} {r.top_patient_share:.0%} "
              f"({r.n_patients_ge20} patients with >= 20 cells); BASE {r.BASE:.0%} PD1 "
              f"{r.PD1:.0%} RTPD1 {r.RTPD1:.0%}; genes {r.median_genes:.0f} (all {med_g:.0f}), "
              f"mt {r.median_pct_mt:.1f}%")
        print(f"  enriched directions (AUROC >= {STATE_DR_AUROC}): {r.n_enriched_directions}")
        for t in long[(long["cluster"] == cl) & (long["rank"] <= STATE_DR_TOP)].itertuples():
            print(f"    {t.direction:7s} AUROC {t.auroc:.2f} (within {t.median_auroc_within:.2f}, "
                  f"{int(t.n_patients_within)}p)  {t.interpretation}")
        print(f"  markers: {r.top_markers}")

    # figures; the cluster x direction heatmap is 03_3's
    # composition of every cluster: patient and treatment
    fig, axes = plt.subplots(1, 2, figsize=(13, 0.34 * len(clusters) + 2.2),
                             gridspec_kw={"width_ratios": [3, 1.2]})
    pats = pd.Series(cohort).value_counts().index.tolist()
    pal = dict(zip(pats, plt.get_cmap("tab20")(np.linspace(0, 1, len(pats)))))
    rows_ = clusters + ["all"]
    for ax, col, order, colours in (
            (axes[0], pd.Series(cohort), pats, pal),
            (axes[1], meta["treatment"].astype(str).reset_index(drop=True),
             ["BASE", "PD1", "RTPD1"], {"BASE": "#2a78d6", "PD1": "#eb6834", "RTPD1": "#1baf7a"})):
        for i, cl in enumerate(rows_):
            m = np.ones(len(lab), bool) if cl == "all" else (lab == cl).values
            share = col[m].value_counts(normalize=True).reindex(order).fillna(0)
            left = 0.0
            for k, v in share.items():
                ax.barh(i, v, left=left, color=colours[k], height=0.7, edgecolor="white",
                        lw=0.4, label=k if i == 0 else None)
                left += v
        ax.set_yticks(range(len(rows_)), [f"{unit} {c}" if c != "all" else "all cells"
                                          for c in rows_], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("share of cells", fontsize=8)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes[0].set_title("patient", fontsize=9)
    axes[1].set_title("treatment", fontsize=9)
    axes[0].legend(fontsize=6, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    axes[1].legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    fig.suptitle(f"Convergence, DRVI states: composition of every {unit}", fontsize=10)
    fig.tight_layout()
    fig.savefig(odir / f"composition_all_clusters_{key}_{drvi_run}.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {odir / f'composition_all_clusters_{key}_{drvi_run}.png'}")

    # Only the partitions with a curated reading.
    if P.get("reading") is not None:
        fig_cluster_factors(summ, long, key, drvi_run, odir, tdir, P["reading"], unit=unit,
                            what=P["what"])


def fig_cluster_factors(summ: pd.DataFrame, long: pd.DataFrame, key: str, drvi_run: str,
                        odir, tdir, reading: dict, unit: str = "cluster",
                        what: str = "DRVI-Leiden cluster") -> None:
    """One bar per cluster: its size (log axis), coloured by its factor, labelled with its DRs."""
    order = list(FACTOR_STYLE)
    df = summ.copy()
    df["factor"] = [reading.get(c, ("none", ""))[0] for c in df.index]
    df["reading"] = [reading.get(c, ("none", "not curated"))[1] for c in df.index]
    df["rank"] = df["factor"].map({f: i for i, f in enumerate(order)})
    df = df.sort_values(["rank", "n_cells"], ascending=[True, False])
    df.drop(columns="rank").to_csv(tdir / f"state_cluster_factors_{key}_{drvi_run}.csv")
    print(f"[table] {tdir / f'state_cluster_factors_{key}_{drvi_run}.csv'}")

    def dr_text(cl):
        top = long[(long["cluster"] == cl) & (long["rank"] <= 2)]
        parts = []
        for r in top.itertuples():
            w = (f", {r.median_auroc_within:.2f} within" if r.n_patients_within
                 and r.n_patients_within > 1 else "")
            parts.append(f"{r.direction} {r.auroc:.2f}{w}")
        return "; ".join(parts)

    n = len(df)
    y = np.arange(n)[::-1]
    fig, ax = plt.subplots(figsize=(13.5, 0.42 * n + 2.4))
    for yi, (cl, r) in zip(y, df.iterrows()):
        colour = FACTOR_STYLE[r["factor"]][0]
        ax.barh(yi, r["n_cells"], height=0.6, color=colour)
        pat = (f"{r['top_patient'].replace('Patient', 'P')} {r['top_patient_share']:.0%}, "
               f"{int(r['n_patients_ge20'])} pts")
        ax.text(r["n_cells"] * 1.12, yi,
                f"{dr_text(cl)}   |   {r['reading']}   |   {pat}",
                va="center", ha="left", fontsize=7, color="0.25")
    ax.set_xscale("log")
    ax.set_xlim(8, df["n_cells"].max() * 400)
    # The right half of the axis holds the labels.
    ticks = [t for t in (10, 100, 1000, 10000, 100000) if t <= df["n_cells"].max() * 2]
    ax.set_xticks(ticks, [f"{t:,}" for t in ticks])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.spines["bottom"].set_bounds(8, ticks[-1])
    ax.set_yticks(y, [f"{unit} {c}" for c in df.index], fontsize=8)
    ax.set_xlabel(f"cells in the {unit} (log scale)", fontsize=9)
    for s_ in ("top", "right", "left"):
        ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_color("0.7")
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="0.9", lw=0.6)
    ax.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=FACTOR_STYLE[f][0]) for f in order
               if f in set(df["factor"])]
    labels = [FACTOR_STYLE[f][1] for f in order if f in set(df["factor"])]
    ax.legend(handles, labels, fontsize=7.5, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.06), ncol=4)
    ax.set_title(f"Convergence, DRVI states: what each {what} is", fontsize=12, loc="left")
    fig.tight_layout()
    fig.savefig(odir / f"cluster_factors_barplot_{key}_{drvi_run}.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {odir / f'cluster_factors_barplot_{key}_{drvi_run}.png'}")


def run_ora(lists: dict[str, list[str]], sets: dict[str, list[str]], background: list[str],
            arm: str, resolution) -> pd.DataFrame:
    """Every gene list against every signature: offline hypergeometric, declared background (BH in apply_bh)."""
    records = []
    for name, genes in lists.items():
        try:
            res = gp.enrich(gene_list=genes, gene_sets=sets, background=background,
                            outdir=None)
        except Exception as exc:                       # noqa: BLE001 - reported, not raised
            print(f"      {name}: ORA FAILED ({exc})", flush=True)
            continue
        res_df = getattr(res, "results", None)
        if res_df is None or len(res_df) == 0:
            continue                                   # no overlap at all: p = 1, added below
        df = pd.DataFrame(res_df).copy()
        df.insert(0, "gene_list", str(name))
        records.append(df)
    if not records:
        return pd.DataFrame(columns=["gene_list", "Term", "P-value", "arm", "resolution"])
    long = pd.concat(records, ignore_index=True)
    long["arm"] = arm
    long["resolution"] = resolution
    return long


def apply_bh(long: pd.DataFrame, n_lists: int, n_sets: int) -> tuple[pd.DataFrame, int]:
    """BH across all n_lists x n_sets pairs of one arm at one resolution; missing pairs added back as p = 1."""
    n_pairs = n_lists * n_sets
    if long.empty:
        return long.assign(fdr_bh=pd.Series(dtype=float),
                           significant=pd.Series(dtype=bool)), n_pairs
    pvals = np.concatenate([long["P-value"].values, np.ones(max(0, n_pairs - len(long)))])
    _, padj, _, _ = multipletests(pvals, alpha=FDR, method="fdr_bh")
    out = long.copy()
    out["fdr_bh"] = padj[:len(long)]
    out["significant"] = out["fdr_bh"] < FDR
    return out, n_pairs




def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smi-resolution", type=float, default=CS.LEIDEN_RESOLUTION,
                   help="the 03_3 Leiden resolution of the clusters, as in leiden_smi_tum.py "
                        "(default %(default)s); one of "
                        f"{list(CS.LEIDEN_RESOLUTIONS)}")
    p.add_argument("--states-de", action="store_true",
                   help="characterise the STATE_PICKS clusters at --smi-resolution and stop: "
                        "one-vs-rest DE on all genes, per-patient consistency of the markers, "
                        "composition, marker overlap with the decoder's top genes, ORA")
    p.add_argument("--states-dr", action="store_true",
                   help="every cluster of --smi-resolution: which live DR directions are "
                        "enriched in it (AUROC, cluster vs rest, and within patient), what each "
                        "direction is (name, Route B hits, 03_3 GO term, confounders, top genes), "
                        "and a one-row summary per cluster. Reads the --states-de DE table")
    p.add_argument("--n-top-genes", type=int, default=N_TOP_GENES,
                   help=f"genes per decoder list for the marker overlap (default {N_TOP_GENES}, "
                        "as 03_7)")
    return p.parse_args()


def main():
    args = parse_args()
    if not (args.states_de or args.states_dr):
        sys.exit("pick a mode: --states-de or --states-dr")
    drvi_run = CS.run_id()
    CS.banner(f"{STEP} - DRVI-Leiden clusters ({drvi_run})")

    if args.states_dr:
        drvi_states_dr(args, drvi_run)
        print("\ndone (--states-dr).")
        return
    hvg = ad.read_h5ad(C.HVG_H5AD)
    hvg.obs_names = hvg.obs_names.astype(str)
    drvi_states_de(hvg, args, drvi_run)
    print("\ndone (--states-de).")


if __name__ == "__main__":
    main()
