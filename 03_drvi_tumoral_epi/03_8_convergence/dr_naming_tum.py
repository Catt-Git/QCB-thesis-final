#!/usr/bin/env python3
"""03_8, the reading step: how many latent dimensions actually get a name.

Reads the three convergence tables together: of the live dimensions (vanished ones pruned),
how many get a name? A name is a convergent row (both routes, same signature family) at
Route A rho >= 0.30 (checked against route_a_null_tum.py's per-direction null) and
Route B FDR < 0.05: a yes or a no, with the rho reported rather than graded. Each dimension's
own correlation with depth, cycle and mito is added, to flag depth axes (DR n is not latent
column n: DRVI titles dimensions by reconstruction_effect).

Figures: programs_named, one bar per program with its max rho plus, in orange, the directions
leiden_smi_tum.py ties to a DRVI-Leiden cluster (STATE_PICKS; --no-states leaves them out);
programs_named_no_rho, the same without rho; dr_names, the same calls read per DR: one row per
dimension, its + and - side, and the name each side gets. A name and a cluster match answer
different questions.

Usage:
    export DATA_DIR=<data root>
    python3 dr_naming_tum.py
    python3 dr_naming_tum.py --all-states
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from types import SimpleNamespace

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402

# The three collections with a convergence table for this run.
COLLECTIONS = ("scie", "emt", "gavish_tnbc")

# The naming bar, the phase-wide Route A bar.
NAMING_RHO_MIN = C.ROUTE_A_RHO_MIN

# A dimension whose own coordinate tracks depth this hard is a depth axis.
DIM_DEPTH_FLAG = 0.30

# Slot 1 of the categorical palette: one category, a named direction.
COL_BAR = "#2a78d6"
# Slot 2, for the second category of the overlay figure: a direction tied to a DRVI state.
COL_STATE = "#eb6834"
# Slot 4 (green), as "sequencing depth" in leiden_states_tum.py.
COL_DEPTH = "#1baf7a"

# Cluster -> direction pairs: the SMI matches at Leiden 1.5, SMI >= 0.7 (all one to one).
# Overlaid here, characterised by 03_8b; re-curate after a re-run of 03_3.
STATE_PICKS = {
    "8": ["DR 6-"],
    "16": ["DR 33+"],
    "20": ["DR 30-"],
    "26": ["DR 11-"],
    "29": ["DR 19-"],
    "30": ["DR 41-"],
    "32": ["DR 44+"],
}

STEP = "03_8_convergence"
# Pseudo-collection: these tables belong to the run, not to a collection.
CONSENSUS = SimpleNamespace(name="consensus",
                            title="all three collections, read together")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--states-resolution", type=float, default=C.CS.LEIDEN_RESOLUTION,
                   help="the DRVI-Leiden resolution whose SMI matches (leiden_smi_tum.py) are "
                        "overlaid on the program bars (default %(default)s)")
    p.add_argument("--all-states", action="store_true",
                   help="overlay every SMI match over the threshold, not only STATE_PICKS")
    p.add_argument("--no-states", action="store_true",
                   help="leave the SMI clusters out of the figures")
    return p.parse_args()


def dim_n(s: str) -> int:
    return int(re.search(r"\d+", str(s)).group())


def read_convergence() -> pd.DataFrame:
    """The three convergence tables, stacked, with the collection as a column."""
    root = C.table_dir(STEP)
    out = []
    for name in COLLECTIONS:
        path = C.table_path("convergence", name)
        if not path.exists():
            print(f"[skip] {path} does not exist")
            continue
        d = pd.read_csv(path, comment="#")
        d["collection"] = name
        out.append(d)
    if not out:
        sys.exit(f"no convergence table under {root}; run 03_8 first")
    d = pd.concat(out, ignore_index=True)
    d["dim_n"] = d["dimension"].map(dim_n)
    return d


def named_rows(conv: pd.DataFrame, rho_log_umi: pd.Series) -> pd.DataFrame:
    """The rows that clear the naming bar, flagged when their dimension is a depth axis."""
    hits = conv[(conv["A_rho"] >= NAMING_RHO_MIN)
                & conv["B_significant"] & conv["same_family"]].copy()
    hits["axis_is_depth"] = hits["dim_n"].map(rho_log_umi).abs() >= DIM_DEPTH_FLAG
    return hits


def dimension_covariates(embed) -> pd.DataFrame:
    """Per dimension: depth, cycle and mito correlation (03_6 does it for the scores), plus DRVI's ranking."""
    X = np.asarray(embed.X)
    obs = embed.obs
    log_umi = np.log10(obs["total_counts"].to_numpy(float) + 1)
    covs = {
        "rho_log_umi": log_umi,
        "rho_n_genes": obs["n_genes_by_counts"].to_numpy(float),
        "rho_S_score": obs["S_score"].to_numpy(float),
        "rho_G2M_score": obs["G2M_score"].to_numpy(float),
        "rho_pct_mt": obs["pct_counts_mt"].to_numpy(float),
    }
    v = embed.var
    rows = []
    for i in range(X.shape[1]):
        r = {"dimension": v["title"].iloc[i],
             "latent_column": int(v["original_dim_id"].iloc[i]),
             "drvi_order": int(v["order"].iloc[i]),
             "reconstruction_effect": float(v["reconstruction_effect"].iloc[i]),
             "vanished": bool(v["vanished"].iloc[i]),
             "vanished_positive_direction": bool(v["vanished_positive_direction"].iloc[i]),
             "vanished_negative_direction": bool(v["vanished_negative_direction"].iloc[i])}
        r.update({k: float(spearmanr(X[:, i], y).statistic) for k, y in covs.items()})
        rows.append(r)
    d = pd.DataFrame(rows)
    d["latent_variance_share"] = d["reconstruction_effect"] / d["reconstruction_effect"].sum()
    d["dim_n"] = d["dimension"].map(dim_n)
    return d.set_index("dim_n").sort_index()


def main():
    args = parse_args()
    C.banner("03_8 - how many dimensions get a name")
    print(f"run {C.RUN_ID}; reading {C.TABLE_DIR.name}/ (vanished tail pruned)")
    print(f"naming bar: Route A rho >= {NAMING_RHO_MIN} AND Route B FDR < 0.05, same family. "
          f"A name is a yes or a no; the magnitude is reported, not graded.")
    print(f"{NAMING_RHO_MIN} against the PER-DIRECTION null of route_a_null_tum.py "
          f"(route_a_null_* tables) - NOT the per-list null of sig_collections")
    print(f"the convergence tables were written at rho >= {C.ROUTE_A_RHO_MIN}, so the bar "
          f"is re-derived here from A_rho / B_significant / same_family: exact, because "
          f"raising a threshold can only remove rows")

    conv = read_convergence()
    adata = ad.read_h5ad(C.EMBED_H5AD)
    cov = dimension_covariates(adata)

    dims = sorted(conv["dim_n"].unique())
    print(f"\n{len(dims)} dimensions x 2 directions = {2 * len(dims)} rows, "
          f"x {conv['collection'].nunique()} collections = {len(conv)} verdicts")

    # Verdict at this bar, rebuilt from the columns convergence_tum.py uses.
    hits = named_rows(conv, cov["rho_log_umi"])
    for c in ("rho_log_umi", "rho_S_score", "rho_G2M_score", "rho_pct_mt",
              "latent_variance_share", "drvi_order", "latent_column"):
        hits[c] = hits["dim_n"].map(cov[c])

    C.banner("the bar")
    stored = conv[conv["verdict"] == "convergent"]
    print(f"convergent at the stored bar ({C.ROUTE_A_RHO_MIN}): {len(stored):3d} rows on "
          f"{stored['dim_n'].nunique():2d} dimensions")
    print(f"convergent at this bar      ({NAMING_RHO_MIN}): {len(hits):3d} rows on "
          f"{hits['dim_n'].nunique():2d} dimensions")

    named = hits.copy()
    clean = named[~named.axis_is_depth]
    print(f"\nnamed:                          {named['dim_n'].nunique()} dimensions, "
          f"{named['dim_direction'].nunique()} directions")
    print(f"  of which the axis is depth:   "
          f"{named[named.axis_is_depth]['dim_n'].nunique()} "
          f"({', '.join(sorted(named[named.axis_is_depth]['dimension'].unique(), key=dim_n))})")
    print(f"named and not a depth axis:     {clean['dim_n'].nunique()} dimensions, "
          f"{clean['B_best_signature'].nunique()} distinct programs")
    two = clean.groupby("dim_direction")["collection"].nunique()
    print(f"  named by two collections:     {(two > 1).sum()} directions "
          f"({', '.join(sorted(two[two > 1].index, key=lambda s: dim_n(s)))})")

    # per-dimension table
    rows = []
    for n in dims:
        sub = hits[hits.dim_n == n].sort_values("A_rho", ascending=False)
        best = sub.iloc[0] if len(sub) else None
        strong = sub
        c = cov.loc[n]
        if len(strong) and abs(c["rho_log_umi"]) >= DIM_DEPTH_FLAG:
            call = "depth_axis"
        elif len(strong):
            call = "named"
        else:
            call = "unnamed"
        rows.append({
            "dimension": f"DR {n}",
            "drvi_order": int(c["drvi_order"]),
            "latent_column": int(c["latent_column"]),
            "latent_variance_share": c["latent_variance_share"],
            "call": call,
            "program": best["B_best_signature"] if best is not None else "",
            "direction": best["dim_direction"] if best is not None else "",
            "A_rho": best["A_rho"] if best is not None else np.nan,
            "B_fdr": best["B_fdr"] if best is not None else np.nan,
            "same_signature": bool(best["same_signature"]) if best is not None else False,
            "named_by": ", ".join(sorted(strong["collection"].unique())),
            "n_convergent_rows": len(sub),
            "other_programs": ", ".join(
                f"{r.dim_direction} {r.B_best_signature} (rho {r.A_rho:.2f})"
                for r in sub.iloc[1:].itertuples()),
            "confounder_flags": best["confounder_flags"] if best is not None else "",
            "rho_log_umi": c["rho_log_umi"], "rho_S_score": c["rho_S_score"],
            "rho_G2M_score": c["rho_G2M_score"], "rho_pct_mt": c["rho_pct_mt"],
        }
        )
    per_dim = pd.DataFrame(rows).set_index("dimension")
    C.write_table(per_dim, "dr_naming", CONSENSUS)

    C.banner("per dimension, DRVI order")
    show = ["drvi_order", "latent_variance_share", "call", "direction", "program",
            "A_rho", "named_by", "rho_log_umi"]
    print(per_dim.sort_values("drvi_order")[show].to_string(float_format="%.3f"))

    # program table
    prog = (clean.groupby("B_best_signature")
            .agg(n_directions=("dim_direction", "nunique"),
                 max_rho=("A_rho", "max"), min_fdr=("B_fdr", "min"),
                 dimensions=("dim_direction",
                             lambda s: ", ".join(sorted(set(s), key=lambda x: dim_n(x)))),
                 collections=("collection", lambda s: ", ".join(sorted(set(s)))))
            .sort_values(["n_directions", "max_rho"], ascending=False))
    prog.index.name = "program"
    C.write_table(prog, "program_dimensions", CONSENSUS)
    C.banner(f"{len(prog)} programs on {clean['dim_n'].nunique()} dimensions")
    print(prog.to_string(float_format="%.3g"))

    # figure
    # programs_named and programs_named_no_rho: the same bars (with the SMI clusters), with and
    # without max rho; dr_names: the same calls read per direction.
    n_named = clean["dim_n"].nunique()
    res = args.states_resolution
    selection = None if args.all_states else STATE_PICKS
    smi = None if args.no_states else read_smi(res, selection)
    if smi is not None:
        smi, prog_of = states_vs_names(smi, clean, named, res, selection is not None)
    for show_rho in (True, False):
        if smi is not None:
            plot_states_overlay(prog, smi, prog_of, n_named, len(dims), res,
                                selection is not None, show_rho=show_rho)
        else:
            plot_program_bars(prog, n_named, len(dims), NAMING_RHO_MIN, show_rho=show_rho)
    plot_dr_names(named, adata, dims, smi, res)


def states_vs_names(smi: pd.DataFrame, clean: pd.DataFrame, named: pd.DataFrame,
                    res: float, selected: bool):
    """What this step says about each SMI-matched direction; writes states_vs_names."""
    prog_of = clean.groupby("dim_direction")["B_best_signature"].agg(lambda s: sorted(set(s)))
    depth_dirs = set(named.loc[named["axis_is_depth"], "dim_direction"])

    def status(d: str) -> str:
        if d in prog_of.index:
            return "named: " + ", ".join(n.replace("_", " ") for n in prog_of[d])
        return "named, depth axis" if d in depth_dirs else "unnamed"

    smi = smi.copy()
    smi["status"] = smi["direction"].map(status)
    smi = smi.sort_values(["cluster", "smi"], ascending=[True, False])
    tbl = smi.rename(columns={"direction": "dim_direction"})
    tbl["resolution"] = res
    C.write_table(tbl, "states_vs_names" + ("" if selected else "_all"), CONSENSUS, index=False)

    state_dirs = set(smi["direction"])
    C.banner(f"DRVI-Leiden cluster matches (resolution {res:g}) against the names")
    print(f"{len(state_dirs)} directions matched to {smi['cluster'].nunique()} clusters; "
          f"{len(state_dirs & set(prog_of.index))} of them also carry a program name")
    print(tbl[["cluster", "dim_direction", "smi", "status"]].to_string(index=False,
                                                                      float_format="%.2f"))
    return smi, prog_of


def plot_program_bars(prog: pd.DataFrame, n_dims_named: int, n_dims: int,
                        rho_min: float, show_rho: bool = True) -> None:
    """One horizontal bar per program, one segment per direction carrying it; max_rho printed
    unless show_rho is False (programs_named_no_rho)."""
    p = prog.iloc[::-1]                       # largest at the top of a horizontal axis
    y = np.arange(len(p))
    fig, ax = plt.subplots(figsize=(10.5, 0.38 * len(p) + 1.9))
    ax.barh(y, p["n_directions"], height=0.46, color=COL_BAR)
    for i, (n, row) in enumerate(p.iterrows()):
        ax.text(row["n_directions"] + 0.10, i,
                row["dimensions"] + (f"   (max rho {row['max_rho']:.2f})" if show_rho else ""),
                va="center", ha="left", fontsize=8, color="0.30")
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ") for n in p.index], fontsize=9)
    ax.set_xlabel("latent directions on which the two routes converge", fontsize=9)
    ax.set_xlim(0, float(p["n_directions"].max()) + 2.6)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.grid(axis="x", color="0.88", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("0.7")
    ax.tick_params(length=0)
    ax.set_title("Convergence: programs named", fontsize=12, loc="left", pad=24)
    ax.text(0, 1.004, f"{n_dims_named} of {n_dims} dimensions named (Route A rho >= {rho_min}, "
            "Route B FDR < 0.05, same family)",
            transform=ax.transAxes, fontsize=8.5, color="0.35", va="bottom")
    fig.tight_layout()
    C.savefig("programs_named" + ("" if show_rho else "_no_rho"), STEP, CONSENSUS, fig=fig)
    plt.close(fig)


def read_smi(res: float, selection: dict[str, list[str]] | None) -> pd.DataFrame | None:
    """leiden_smi_tum.py's SMI matches, cut to `selection`; None if the table is missing."""
    path = C.table_dir(STEP) / f"smi_matches_drvi_leiden_{res:g}_{C.RUN_ID}.csv"
    if not path.exists():
        print(f"\n[skip] no SMI table at {path}: run leiden_smi_tum.py --smi-resolution {res:g}")
        return None
    smi = pd.read_csv(path)
    smi["cluster"] = smi["cluster"].astype(str)
    if selection is not None:
        keep = {(c, d) for c, ds in selection.items() for d in ds}
        missing = keep - set(zip(smi["cluster"], smi["direction"]))
        if missing:
            sys.exit(f"selected pairs not in {path.name}: {sorted(missing)}")
        smi = smi[[(c, d) in keep for c, d in zip(smi["cluster"], smi["direction"])]]
    return smi


def plot_states_overlay(prog: pd.DataFrame, smi: pd.DataFrame, prog_of: pd.Series,
                        n_dims_named: int, n_dims: int, res: float, selected: bool,
                        show_rho: bool = True) -> None:
    """The program bars plus the directions leiden_smi_tum.py matched to a DRVI-Leiden cluster."""
    state_dirs = set(smi["direction"])
    cl_of = smi.groupby("direction")["cluster"].agg(lambda s: ", ".join(sorted(s, key=int)))

    # ---- geometry: programs on top, a gap, then one bar per cluster
    p = prog.copy()
    p["n_state"] = p["dimensions"].map(
        lambda s: sum(d.strip() in state_dirs for d in s.split(",")))
    clusters = (smi.groupby("cluster")
                .agg(n=("direction", "nunique"), top=("smi", "max"),
                     label=("direction", lambda s: ", ".join(s)))
                .sort_values(["n", "top"], ascending=False))
    n_p, n_c = len(p), len(clusters)
    n_shared = len(state_dirs & set(prog_of.index))
    y_p = np.arange(n_p)[::-1] + n_c + 1.2           # programs, largest at the top
    y_c = np.arange(n_c)[::-1]                        # clusters below
    xmax = float(max(p["n_directions"].max(), clusters["n"].max())) + 3.2

    fig, ax = plt.subplots(figsize=(11.5, 0.36 * (n_p + n_c) + 2.6))
    ax.barh(y_p, p["n_directions"], height=0.46, color=COL_BAR,
            label="named direction (both routes, this step)")
    ax.barh(y_p, p["n_state"], height=0.46, color=COL_STATE,
            label=f"direction matched to a DRVI-Leiden cluster (SMI, resolution {res:g}"
                  + (", selected pairs)" if selected else ", all >= 0.7)"))
    for yi, (_, row) in zip(y_p, p.iterrows()):
        parts = [f"{d.strip()} [cl {cl_of[d.strip()]}]" if d.strip() in state_dirs else d.strip()
                 for d in row["dimensions"].split(",")]
        ax.text(row["n_directions"] + 0.10, yi,
                ", ".join(parts) + (f"   (max rho {row['max_rho']:.2f})" if show_rho else ""),
                va="center", ha="left", fontsize=7.5, color="0.30")
    ax.barh(y_c, clusters["n"], height=0.46, color=COL_STATE)
    for yi, (cl, row) in zip(y_c, clusters.iterrows()):
        sub = smi[smi["cluster"] == cl]
        # The SMI score goes with rho: the _no_rho figure shows only the direction.
        def label(r) -> str:
            parts = ([f"{r.smi:.2f}"] if show_rho else []) + \
                ([] if r.status == "unnamed" else [r.status])
            return r.direction + (f" ({'; '.join(parts)})" if parts else "")

        txt = ", ".join(label(r) for r in sub.itertuples())
        ax.text(row["n"] + 0.10, yi, txt, va="center", ha="left", fontsize=7.5, color="0.30")

    ax.axhline(n_c + 0.1, color="0.7", lw=0.8)
    ax.text(-0.02, y_p[0] + 0.75, "programs named by both routes", transform=ax.get_yaxis_transform(),
            ha="right", fontsize=8.5, color="0.35", style="italic")
    ax.text(-0.02, y_c[0] + 0.75, f"DRVI-Leiden clusters ({res:g})",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=8.5, color="0.35",
            style="italic")
    ax.set_yticks(np.concatenate([y_p, y_c]))
    ax.set_yticklabels([n.replace("_", " ") for n in p.index]
                       + [f"cluster {c}" for c in clusters.index], fontsize=8.5)
    ax.set_xlabel("latent directions", fontsize=9)
    ax.set_xlim(0, xmax)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.grid(axis="x", color="0.88", lw=0.6)
    ax.set_axisbelow(True)
    for s_ in ("top", "right", "left"):
        ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_color("0.7")
    ax.tick_params(length=0)
    ax.legend(loc="lower right", fontsize=7.5, frameon=False)
    ax.set_title("Convergence: programs named, with the DRVI states",
                 fontsize=12, loc="left", pad=24)
    ax.text(0, 1.004, f"blue = named ({n_dims_named} of {n_dims} dimensions); orange = the "
            f"direction of a DRVI-Leiden cluster (SMI), {n_shared} of {len(state_dirs)} named",
            transform=ax.transAxes, fontsize=8.2, color="0.35", va="bottom")
    fig.tight_layout()
    C.savefig("programs_named" + ("" if show_rho else "_no_rho")
              + ("" if selected else "_all"), STEP, CONSENSUS, fig=fig)
    plt.close(fig)


def plot_dr_names(named: pd.DataFrame, embed, dims: list[int], smi: pd.DataFrame | None,
                  res: float) -> None:
    """One row per live dimension, its + and - side as two tiles: the name each side gets.

    With `smi` a direction matched to a DRVI-Leiden cluster carries its cluster too.
    """
    var = embed.var.set_index(embed.var["title"].astype(str))
    progs = named.groupby("dim_direction")["B_best_signature"].agg(lambda s: sorted(set(s)))
    depth = set(named.loc[named["axis_is_depth"], "dim_direction"])
    cluster = ({} if smi is None else
               smi.groupby("direction")["cluster"].agg(lambda s: ", ".join(sorted(s, key=int))))
    vanished = {"+": "vanished_positive_direction", "-": "vanished_negative_direction"}
    colour = {"named": COL_BAR, "depth": COL_DEPTH, "cluster": COL_STATE,
              "unnamed": "#ececec", "vanished": "white"}

    def call(d: str, side: str) -> tuple[str, str]:
        key = f"DR {d}{side}"
        if bool(var.loc[f"DR {d}", vanished[side]]):
            return "vanished", "vanished"
        cl = f"  [cluster {cluster[key]}]" if key in cluster else ""
        if key in progs.index:
            txt = " / ".join(n.replace("_", " ") for n in progs[key])
            return ("depth" if key in depth else "named"), txt + cl
        if cl:
            return "cluster", cl.strip()[1:-1]
        return "unnamed", ""

    fig, ax = plt.subplots(figsize=(11.5, 0.24 * len(dims) + 1.9))
    counts = {k: 0 for k in colour}
    for i, d in enumerate(dims):
        y = len(dims) - 1 - i
        for j, side in enumerate(("+", "-")):
            kind, txt = call(d, side)
            counts[kind] += 1
            ax.add_patch(plt.Rectangle((j, y - 0.42), 0.97, 0.84, color=colour[kind],
                                       ec="0.85" if kind == "vanished" else colour[kind], lw=0.6))
            ink = "white" if kind in ("named", "depth", "cluster") else "0.55"
            ax.text(j + 0.02, y, txt, va="center", ha="left", fontsize=6.3, color=ink)
    ax.set_xlim(0, 1.97)
    ax.set_ylim(-0.7, len(dims) - 0.3)
    ax.set_yticks(range(len(dims)), [f"DR {d}" for d in dims][::-1], fontsize=7)
    ax.set_xticks([0.485, 1.485], ["+ side", "- side"], fontsize=9)
    ax.xaxis.tick_top()
    for s_ in ("top", "right", "left", "bottom"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(length=0)
    n_live = 2 * len(dims) - counts["vanished"]
    labels = {"named": "named (both routes)",
              "depth": "named, but the dimension is a depth axis",
              "cluster": f"DRVI-Leiden cluster only (SMI, {res:g})" if smi is not None else "",
              "unnamed": "no name"}
    keys = [k for k in labels if labels[k]]
    handles = [plt.Rectangle((0, 0), 1, 1, color=colour[k]) for k in keys]
    ax.legend(handles, [labels[k] for k in keys],
              loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=len(keys), fontsize=7.5,
              frameon=False)
    # Title and subtitle in figure coordinates: the side labels sit on top of the axes.
    fig.text(0.01, 0.995, "Convergence: the name of each DR", fontsize=12, va="top")
    sub = (f"{counts['named'] + counts['depth']} of {n_live} live directions "
           f"named ({counts['depth']} on a depth axis)")
    if smi is not None:
        sub += f", {counts['cluster']} more tied to a cluster only"
    fig.text(0.01, 0.975, sub, fontsize=8.5, color="0.35", va="top")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    C.savefig("dr_names", STEP, CONSENSUS, fig=fig)
    plt.close(fig)


if __name__ == "__main__":
    main()
