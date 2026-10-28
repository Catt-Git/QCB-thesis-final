#!/usr/bin/env python3
"""03_8, the picture of a name: the program and its axis, on the same UMAP.

One figure per bar segment of dr_naming_tum.py's programs_named (pairs read from
program_dimensions_consensus_*.csv): left, the within-cohort z of the program score
(--score raw for the raw one); right, the coordinate on the dimension times the sign of the
direction. Both one-sided: grey at or below zero, a ramp above (clipped at p99), red for a +
and blue for a - direction. If the name holds, the coloured cells match. UMAP of the DRVI
latent (03_3). Rho and FDR are read from the 03_8 tables; only the Spearman of the two plotted
vectors is computed here.

Usage:
    export DATA_DIR=<data root>
    python3 dr_umap_tum.py
    python3 dr_umap_tum.py --only "DR 2+,DR 15-,DR 12-"
    python3 dr_umap_tum.py --score raw
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
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C          # noqa: E402
import sig_collections as SC          # noqa: E402

# The naming bar, imported from the step that applies it.
from dr_naming_tum import NAMING_RHO_MIN as NAMING_BAR   # noqa: E402

# One figure per bar: a subfolder, so they do not bury the step's other figures.
STEP = "03_8_convergence/dr_umap"

# Same pseudo-collection as dr_naming.
CONSENSUS = SimpleNamespace(name="consensus",
                            title="all three collections, read together")

# Flat grey at or below zero, one saturated ramp above it.
GREY = "#dcdcdf"

RAMPS = {
    "score": (GREY, "#b88cf0", "#8b2fe0", "#5b00b0"),   # the program: purple
    "+":     (GREY, "#f59b86", "#e5352a", "#a00a00"),   # a positive direction: red
    "-":     (GREY, "#8fc0ef", "#1f6fd0", "#0a3577"),   # a negative direction: blue
}

# Colour scale clipped at this quantile of the positive values only.
CLIP_Q = 0.99

POINT_SIZE = 2.4


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--only", default="",
                   help="comma-separated directions to draw, e.g. 'DR 50+,DR 4+,DR 9+' "
                        "(default: every pair in the barplot)")
    p.add_argument("--program", default="",
                   help="comma-separated program names to draw, e.g. 'MP7_STRESS_IN_VITRO'")
    p.add_argument("--score", choices=("z", "raw"), default="z",
                   help="the left panel: within-cohort z (default, what Route A correlates) "
                        "or the raw score_genes value")
    p.add_argument("--dpi", type=int, default=200,
                   help="200 by default and not the phase's 300: these are 41,943-point "
                        "scatters and there are 32 of them")
    return p.parse_args()


def dim_n(s: str) -> int:
    return int(re.search(r"\d+", str(s)).group())


def read_pairs(args) -> pd.DataFrame:
    """The bars of programs_named, one row per segment (program x direction x collection)."""
    prog = C.read_table("program_dimensions", CONSENSUS)
    rows = []
    for program, r in prog.iterrows():
        for direction in [d.strip() for d in str(r["dimensions"]).split(",")]:
            for coll in [c.strip() for c in str(r["collections"]).split(",")]:
                rows.append({"program": program, "dim_direction": direction,
                             "collection": coll, "dimension": direction[:-1].strip(),
                             "sign": 1.0 if direction.endswith("+") else -1.0,
                             "max_rho": float(r["max_rho"])})
    d = pd.DataFrame(rows)

    if args.only:
        want = {s.strip() for s in args.only.split(",")}
        unknown = want - set(d["dim_direction"])
        if unknown:
            print(f"[note] not on the barplot, ignored: {', '.join(sorted(unknown))}")
        d = d[d["dim_direction"].isin(want)]
    if args.program:
        want = {s.strip() for s in args.program.split(",")}
        d = d[d["program"].isin(want)]
    if d.empty:
        sys.exit("no pair left to draw after --only/--program")
    return d.sort_values(["dim_direction", "program"],
                         key=lambda s: s.map(C.dim_sort_key) if s.name == "dim_direction" else s)


def collection_by_slug(slug: str):
    """A collection by the slug its outputs are written under (gavish_tnbc is not a --collection value)."""
    for c in (*SC.COLLECTIONS.values(), *SC.GAVISH_VARIANTS.values()):
        if c.name == slug:
            return c
    raise KeyError(f"no collection writes under the slug {slug!r}")


def load_collection(name: str, which: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per-cell scores, the convergence table and the Route A rho matrix, for one collection."""
    coll = collection_by_slug(name)
    path = C.scores_csv(coll)
    if not path.exists():
        sys.exit(f"missing {path}: run 03_6_cell_first/cell_first_tum.py --collection {name}")
    scores = pd.read_csv(path, index_col=0, comment="#")
    prefix = "z_" if which == "z" else "score_"
    keep = [c for c in scores.columns if c.startswith(prefix)]
    scores = scores[keep]
    scores.columns = [c[len(prefix):] for c in scores.columns]

    conv = C.read_table("convergence", coll)
    rho = C.read_table("dim_signature_spearman", coll)
    print(f"[read] {name:12s} {scores.shape[0]} cells x {scores.shape[1]} readouts "
          f"({prefix.rstrip('_')}), {len(conv)} convergence rows")
    return scores, conv, rho


def dimension_frame(adata) -> pd.DataFrame:
    """Cells x dimensions, columns titled 'DR 1' ... from var['title'] (DR n is not column n)."""
    return pd.DataFrame(np.asarray(adata.X), index=adata.obs_names,
                        columns=adata.var["title"].values)


def top_limit(v: np.ndarray) -> float:
    """The high end of the colour scale: the 99th percentile of the POSITIVE values."""
    pos = v[v > 0]
    if pos.size == 0:
        return 1.0
    return float(np.quantile(pos, CLIP_Q)) or float(pos.max()) or 1.0


def panel(ax, umap: np.ndarray, v: np.ndarray, ramp: str, label: str, tick_sign: float = 1.0):
    """One UMAP, grey below zero, `ramp` above it, the strongest cells drawn last.

    `tick_sign` = -1 labels the colourbar of a negated coordinate (negative direction).
    """
    cmap = LinearSegmentedColormap.from_list(f"one_sided_{ramp}", RAMPS[ramp])
    order = np.argsort(v)
    s = ax.scatter(umap[order, 0], umap[order, 1], c=v[order], cmap=cmap,
                   vmin=0.0, vmax=top_limit(v), s=POINT_SIZE, lw=0, rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("0.85")
    ax.set_xlabel("UMAP 1", fontsize=8, color="0.35")
    ax.set_ylabel("UMAP 2", fontsize=8, color="0.35")
    cb = ax.figure.colorbar(s, ax=ax, fraction=0.046, pad=0.02, extend="max")
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    if tick_sign < 0:
        cb.ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda t, _: f"{-t if t else 0.0:.1f}"))
    cb.outline.set_visible(False)
    return s


def draw_pair(row, umap, L, scores, conv, rho_mat, args) -> None:
    program, direction = row["program"], row["dim_direction"]
    dim, sign = row["dimension"], row["sign"]

    if program not in scores.columns:
        print(f"[skip] {direction} {program}: no {args.score} column in "
              f"{row['collection']}'s score table")
        return
    # Clipped at zero; the coordinate is first multiplied by the sign of the direction.
    y = np.clip(scores[program].to_numpy(float), 0.0, None)
    x = np.clip((sign * L[dim]).to_numpy(float), 0.0, None)

    # Read, not measured: the two routes' own numbers for this direction.
    c = conv.loc[direction]
    # Route A rho of the plotted signature, oriented; A_rho may belong to a sibling signature.
    rho_plotted = float(sign * rho_mat.loc[dim, program]) if program in rho_mat.columns \
        else float("nan")
    # Recomputed on the plotted vectors, as a check on the caption.
    rho_drawn = float(spearmanr(sign * L[dim].to_numpy(float),
                                scores[program].to_numpy(float)).statistic)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0))

    score_label = ("within-cohort z" if args.score == "z" else "score_genes") \
        + f", {program.replace('_', ' ')}   (<= 0 grey)"
    panel(axes[0], umap, y, "score", score_label)
    axes[0].set_title(f"cells: {program.replace('_', ' ')}",
                      fontsize=10, loc="left", color="0.2")

    side = "+" if sign > 0 else "-"
    axis_label = (f"DRVI coordinate on {dim}   "
                  f"({'<= 0' if sign > 0 else '>= 0'} grey)")
    panel(axes[1], umap, x, side, axis_label, tick_sign=sign)
    axes[1].set_title(f"axis: {direction}", fontsize=10, loc="left", color="0.2")

    flags = str(c["confounder_flags"])
    bits = [f"Route A rho {rho_plotted:+.3f} on this signature",
            f"Route B FDR {float(c['B_fdr']):.1e}",
            f"named out of {row['collection']}"]
    if str(c["A_best_signature"]) != program:
        bits.append(f"Route A's own best here: {c['A_best_signature']} "
                    f"(rho {float(c['A_rho']):+.3f})")
    if flags and flags != "none":
        bits.append(f"flags: {flags}")

    fig.suptitle(f"Convergence: {program.replace('_', ' ')} on {direction}", fontsize=13,
                 x=0.012, ha="left", y=1.02)
    fig.text(0.012, 0.985, " · ".join(bits), fontsize=8, color="0.35", ha="left", va="top")
    print(f"  {direction} {program}: Spearman of the plotted vectors {rho_drawn:+.3f}")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    C.savefig(f"umap_{C.dim_slug(direction)}_{program}", STEP, CONSENSUS, fig=fig,
              dpi=args.dpi)
    plt.close(fig)

    return {"dim_direction": direction, "program": program,
            "collection": row["collection"], "rho_this_signature": rho_plotted,
            "A_best_signature": str(c["A_best_signature"]), "A_rho": float(c["A_rho"]),
            "B_fdr": float(c["B_fdr"]), "confounder_flags": flags}


def main():
    args = parse_args()
    C.banner("03_8 - the program and its axis on the same UMAP")
    print(f"run {C.RUN_ID}; reading {C.TABLE_DIR.name}/ (vanished tail pruned)")

    pairs = read_pairs(args)
    print(f"{len(pairs)} program-direction pairs on "
          f"{pairs['dim_direction'].nunique()} directions "
          f"({pairs['dimension'].nunique()} dimensions), "
          f"from {', '.join(sorted(pairs['collection'].unique()))}")

    adata = ad.read_h5ad(C.EMBED_H5AD)
    if "X_umap" not in adata.obsm:
        sys.exit(f"{C.EMBED_H5AD.name} has no X_umap: it was written before 03_3's "
                 f"neighbors/umap step")
    umap = np.asarray(adata.obsm["X_umap"])
    L = dimension_frame(adata)
    print(f"[read] {C.EMBED_H5AD.name}: {adata.n_obs} cells, {adata.n_vars} dimensions, "
          f"UMAP on the latent itself")

    loaded = {}
    for name in sorted(pairs["collection"].unique()):
        scores, conv, rho_mat = load_collection(name, args.score)
        missing = L.index.difference(scores.index)
        if len(missing):
            sys.exit(f"{name}: {len(missing)} cells of the embedding have no score row")
        loaded[name] = (scores.loc[L.index], conv, rho_mat)

    C.banner(f"{len(pairs)} figures")
    drawn = []
    for _, row in pairs.iterrows():
        scores, conv, rho_mat = loaded[row["collection"]]
        rec = draw_pair(row, umap, L, scores, conv, rho_mat, args)
        if rec is not None:
            drawn.append(rec)

    report(pd.DataFrame(drawn))
    print(f"\n[done] {C.fig_dir(STEP)}")


def report(d: pd.DataFrame) -> None:
    """Print the rho of the program on each bar, which can differ from max_rho (a sibling signature)."""
    if d.empty:
        return
    d = d.copy()
    d["same"] = d["A_best_signature"] == d["program"]
    C.banner("the program on the bar, against the axis")
    print(d.sort_values("rho_this_signature", ascending=False)
          [["dim_direction", "program", "collection", "rho_this_signature",
            "A_best_signature", "A_rho", "confounder_flags"]]
          .to_string(index=False, float_format="%.3f"))
    split = d[~d["same"]]
    print(f"\n{len(split)} of {len(d)} bars are labelled with a signature that is not the "
          f"one Route A scored highest on that direction")
    weak = split[split["rho_this_signature"] < NAMING_BAR]
    if len(weak):
        print(f"of those, {len(weak)} where the labelled program's OWN rho is BELOW the "
              f"naming bar ({NAMING_BAR}), i.e. the direction is over the bar on a sibling "
              f"signature of the same family and not on the one the bar is labelled with:")
        for r in weak.sort_values("rho_this_signature").itertuples():
            print(f"  {r.dim_direction:8s} {r.program:28s} rho {r.rho_this_signature:+.3f}"
                  f"   (bar quotes {r.A_best_signature} at {r.A_rho:+.3f})")


if __name__ == "__main__":
    main()
