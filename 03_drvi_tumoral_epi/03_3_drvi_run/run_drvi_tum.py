#!/usr/bin/env python
"""
03_3 run_drvi_tum: DRVI on the malignant cells, headless (the compute half of drvi_tum_64.ipynb)

DRVI with 64 latent dimensions on the 2,000 HVGs without MT-, batch_key='cohort',
400 epochs without early stopping, no labels_key. The cells are then clustered with Leiden
on the live DRVI dimensions (k = 15, resolutions 0.25 - 2.0; 1.5 is drvi_leiden); clusters
under 100 cells are merged into their neighbours.

Input : $DATA_DIR/03_tum/shiao_tum_hvg_2k_nomt.h5ad
        Output of 03_2 (reduce_data_tum.py); .layers['counts'] = raw counts.
        $DATA_DIR/03_tum/shiao_tum.h5ad
        Definitive object of 03_2, for the downstream object.
        
Output: $DATA_DIR/03_tum/model_drvi_tum_64_nomt.pt
        $DATA_DIR/03_tum/embed_drvi_tum_64_nomt.h5ad
        Latent space, per-dimension stats, OOD/IND interpretability scores, the Leiden
        columns drvi_leiden_<res> and drvi_leiden (read by 03_4 on).
        $DATA_DIR/03_tum/shiao_tum_drvi_tum_64_nomt.h5ad
        shiao_tum.h5ad + obsm['X_drvi'] + the Leiden columns, for 03_5 and 03_6.
        ../tables/03_3_drvi_tum_64_nomt/drvi_leiden_clusters_drvi_tum_64_nomt.csv
        Clusters per resolution, before and after merging the small ones.
        ../figures/03_3_drvi_tum_64_nomt/
        The notebook's figures; a plotting failure does not fail the run.

Existing outputs are reused; --overwrite recomputes everything.

Usage:
  export DATA_DIR=<data root>
  python3 run_drvi_tum.py
  python3 run_drvi_tum.py --overwrite      # retrain and rewrite everything
  python3 run_drvi_tum.py --early-stopping # stop early instead (see --epochs)
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scvi  # noqa: E402
import drvi  # noqa: E402
from drvi.model import DRVI  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import cell_set as C  # noqa: E402

BATCH_KEY = C.BATCH_KEY
LABEL_KEY = C.LABEL_KEY          # constant 'malignant' here
LEIDEN_KEY = C.LEIDEN_KEY

# The notebook's UMAP panels; constant columns are skipped at draw time.
UMAP_QC_KEYS = {
    "cell_type": LABEL_KEY,
    "leiden": LEIDEN_KEY,
    "cohort": BATCH_KEY,
    "treatment": "treatment",
    "response": "response",
    "phase": "phase",
    "cnv_score": "cnv_score",
    "cnv_corr": "cnv_corr",
    "n_genes_by_counts": "n_genes_by_counts",
    "total_counts": "total_counts",
    "mito": "pct_counts_mt",
    "ribo": "pct_counts_ribo",
    "size_factors": "size_factors",
}
UMAP_COMBINED_KEYS = [LABEL_KEY, LEIDEN_KEY, BATCH_KEY, "treatment", "response", "phase"]

# Groupings of the heatmaps, each drawn sorted by category.
HEATMAP_KEYS = [LEIDEN_KEY, BATCH_KEY, "treatment", "response", "phase", LABEL_KEY]

# Groupings that get one UMAP panel per level.
PANEL_KEYS = [LEIDEN_KEY]

# Every random step of 03_3: training, kNN graph, UMAP, Leiden, plot order.
SEED = 123
UMAP_SEED = SEED



def parse_args():
    p = argparse.ArgumentParser(description="Run DRVI on the malignant epithelial subset")
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--epochs", type=int, default=400,
                   help="epochs; the KL warmup spans all of them [default: 400]")
    p.add_argument("--early-stopping", action="store_true",
                   help="stop after 50 epochs without improvement. Off by default: a run "
                        "that stops before the end of the KL warmup never reaches full KL "
                        "weight, which is what makes unused dimensions vanish")
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR"),
                   help="directory holding the datasets [default: $DATA_DIR]")
    p.add_argument("--fig-dir", default=None,
                   help="where the figures go [default: the repo's "
                        "03_drvi_tumoral_epi/figures/03_3_drvi_tum_64_nomt]")
    p.add_argument("--overwrite", action="store_true",
                   help="retrain and rewrite everything instead of reusing what is on disk")
    return p.parse_args()


def save_model(model, path):
    """Save a DRVI model as the single file `path` (model.save() only writes <dir>/model.pt)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent) as scratch:
        model.save(scratch, overwrite=True)
        Path(scratch, "model.pt").replace(path)
    return path


def load_model(path, adata):
    """Load a model saved by save_model, via a scratch dir holding a symlink to it."""
    path = Path(path).resolve()
    with tempfile.TemporaryDirectory(dir=path.parent) as scratch:
        os.symlink(path, Path(scratch, "model.pt"))
        return DRVI.load(scratch, adata)


def savefig(fig_dir, run_id, name, fig=None, dpi=300):
    """Save and close a figure (the current one by default) as <name>_<run_id>.png.

    For the DRVI and seaborn plots, which ignore sc.settings.figdir.
    """
    fig = plt.gcf() if fig is None else fig
    path = Path(fig_dir) / f"{name}_{run_id}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path}", flush=True)
    return path


def merge_small_clusters(lab: np.ndarray, conn, min_cells: int = C.LEIDEN_MIN_CELLS) -> np.ndarray:
    """Cells of clusters under `min_cells` go to the cluster that holds most of their kNN weight.

    Repeated until no small cluster has a neighbour outside the small clusters; a small cluster
    with no edge to the rest of the graph is kept as it is. Labels are renumbered by size.
    """
    lab = lab.copy()
    conn = conn.tocsr()
    while True:
        sizes = pd.Series(lab).value_counts()
        small = sizes.index[sizes < min_cells]
        keep = ~np.isin(lab, small)
        moved = 0
        for i in np.flatnonzero(~keep):
            nb = conn.indices[conn.indptr[i]:conn.indptr[i + 1]]
            w = conn.data[conn.indptr[i]:conn.indptr[i + 1]]
            ok = keep[nb]
            if ok.any():
                lab[i] = pd.Series(w[ok]).groupby(lab[nb[ok]]).sum().idxmax()
                moved += 1
        if moved == 0:
            break
    order = pd.Series(lab).value_counts().index
    return pd.Series(lab).map({c: k for k, c in enumerate(order)}).to_numpy()


def leiden_params() -> dict:
    """The settings the stored Leiden columns were made with (kept in embed.uns)."""
    return {"k": C.LEIDEN_K, "resolutions": [float(r) for r in C.LEIDEN_RESOLUTIONS],
            "resolution": float(C.LEIDEN_RESOLUTION), "min_cells": C.LEIDEN_MIN_CELLS,
            "seed": SEED}


def leiden_is_current(embed) -> bool:
    """True when the embedding holds every Leiden column, made with the current settings."""
    keys = [C.leiden_key(r) for r in C.LEIDEN_RESOLUTIONS] + [LEIDEN_KEY]
    stored = embed.uns.get("drvi_leiden_params")
    if stored is None or not all(k in embed.obs for k in keys):
        return False
    stored = {k: (list(map(float, v)) if k == "resolutions" else v) for k, v in stored.items()}
    return stored == leiden_params()


def live_umap(embed):
    """Neighbours and UMAP on the live dimensions only (obsm['X_live']), the graph of the Leiden."""
    live = ~embed.var["vanished"].astype(bool).to_numpy()
    embed.obsm["X_live"] = np.asarray(embed.X)[:, live].astype(np.float32)
    sc.pp.neighbors(embed, n_neighbors=C.LEIDEN_K, use_rep="X_live", random_state=SEED)
    sc.tl.umap(embed, random_state=SEED)
    print(f"UMAP on the {int(live.sum())} live dimensions", flush=True)


def leiden_graph(embed, seed=SEED) -> ad.AnnData:
    """kNN graph (k = C.LEIDEN_K) on the live DRVI dimensions, the one every Leiden run uses."""
    live = ~embed.var["vanished"].astype(bool).to_numpy()
    g = ad.AnnData(obs=pd.DataFrame(index=embed.obs_names))
    g.obsm["X_live"] = np.asarray(embed.X)[:, live].astype(np.float32)
    sc.pp.neighbors(g, n_neighbors=C.LEIDEN_K, use_rep="X_live", random_state=seed)
    print(f"Leiden graph: k = {C.LEIDEN_K} on {int(live.sum())} live dimensions", flush=True)
    return g


def leiden_labels(g, res, seed=SEED) -> tuple[np.ndarray, np.ndarray]:
    """Leiden at `res`, then merge_small_clusters: (raw labels, merged labels)."""
    sc.tl.leiden(g, resolution=res, key_added="_leiden", flavor="igraph", n_iterations=2,
                 directed=False, random_state=seed)
    raw = g.obs["_leiden"].astype(int).to_numpy()
    return raw, merge_small_clusters(raw, g.obsp["connectivities"])


def drvi_leiden(embed, seed=SEED) -> pd.DataFrame:
    """Leiden on the live DRVI dimensions at every resolution of C.LEIDEN_RESOLUTIONS.

    Clusters under C.LEIDEN_MIN_CELLS are merged into their neighbours (merge_small_clusters).
    Adds obs drvi_leiden_<res> and drvi_leiden; returns the cluster sizes per resolution.
    """
    g = leiden_graph(embed, seed)
    rows = []
    for res in C.LEIDEN_RESOLUTIONS:
        raw, lab = leiden_labels(g, res, seed)
        embed.obs[C.leiden_key(res)] = pd.Categorical(
            lab.astype(str), categories=[str(c) for c in sorted(np.unique(lab))])
        n = pd.Series(lab).value_counts()
        raw_sizes = pd.Series(raw).value_counts()
        rows.append({"resolution": res, "n_clusters_leiden": len(raw_sizes), "n_clusters": len(n),
                     "cells_in_small_clusters":
                         int(raw_sizes[raw_sizes < C.LEIDEN_MIN_CELLS].sum()),
                     "min_cells": int(n.min()), "median_cells": float(n.median()),
                     "max_cells": int(n.max())})
        print(f"  resolution {res:<5} -> {len(raw_sizes):3d} clusters, {len(n):3d} after merging "
              f"those under {C.LEIDEN_MIN_CELLS} (smallest {n.min():,}, largest {n.max():,})",
              flush=True)
    embed.obs[LEIDEN_KEY] = embed.obs[C.leiden_key(C.LEIDEN_RESOLUTION)]
    embed.uns["drvi_leiden_params"] = leiden_params()
    return pd.DataFrame(rows)


def leiden_resolution_grid(embed, fig_dir, run_id):
    """Every Leiden resolution on the DRVI UMAP, the one used marked."""
    n = len(C.LEIDEN_RESOLUTIONS)
    ncol = 4 if n > 4 else n
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.4 * nrow), squeeze=False)
    for ax, res in zip(axes.flat, C.LEIDEN_RESOLUTIONS):
        key = C.leiden_key(res)
        sc.pl.umap(embed, color=key, ax=ax, show=False, frameon=False, legend_loc="on data",
                   legend_fontsize=5, legend_fontoutline=1.5, size=4,
                   title=f"resolution {res} - {embed.obs[key].nunique()} clusters"
                         + (" (used)" if res == C.LEIDEN_RESOLUTION else ""))
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    n_live = int((~embed.var["vanished"].astype(bool)).sum())
    fig.suptitle(f"DRVI {run_id}: Leiden on the {n_live} live dimensions, on their UMAP",
                 fontsize=11)
    savefig(fig_dir, run_id, "umap_drvi_leiden_all_resolutions", fig=fig)


def groupable(embed, keys):
    """The keys present in embed.obs with more than one level (cell_type is constant here)."""
    out = []
    for key in keys:
        if key not in embed.obs:
            print(f"[skip] grouping {key!r}: not an obs column", flush=True)
        elif embed.obs[key].astype(str).nunique() < 2:
            print(f"[skip] grouping {key!r}: constant "
                  f"({embed.obs[key].astype(str).iloc[0]!r})", flush=True)
        else:
            out.append(key)
    return out


def draw_figures(model, adata, embed, fig_dir, run_id, heatmap_keys):
    """The plotting section of the notebook, in order and with the same names."""
    sc.settings.figdir = fig_dir

    # Same seeded permutation as the notebook, so overplotting hides the same cells.
    order = np.random.default_rng(UMAP_SEED).permutation(embed.n_obs)
    embed_plot = ad.AnnData(
        obs=embed.obs.iloc[order].copy(),
        obsm={"X_umap": embed.obsm["X_umap"][order]},
        # Same palettes as the 03_2 UMAPs.
        uns={k: v for k, v in adata.uns.items() if k.endswith("_colors")},
    )

    for label, col in UMAP_QC_KEYS.items():
        if col not in embed_plot.obs:
            print(f"[skip] umap {label}: no obs column {col!r}", flush=True)
            continue
        if embed_plot.obs[col].dtype.kind not in "fiu" and \
                embed_plot.obs[col].astype(str).nunique() < 2:
            print(f"[skip] umap {label}: {col!r} is constant "
                  f"({embed_plot.obs[col].astype(str).iloc[0]!r})", flush=True)
            continue
        sc.pl.umap(embed_plot, color=col, show=False, save=f"_{label}_{run_id}.png")
        plt.close("all")

    combined = [k for k in UMAP_COMBINED_KEYS
                if k in embed_plot.obs and embed_plot.obs[k].astype(str).nunique() > 1]
    with plt.rc_context({"figure.figsize": (7, 7)}):
        sc.pl.umap(embed_plot, color=combined, ncols=2, wspace=0.8,
                   hspace=0.25, show=False, save=f"_combined_{run_id}.png")
        plt.close("all")

    # One panel per level of each grouping, each highlighting only its own cells.
    for key in PANEL_KEYS:
        if key not in embed_plot.obs or embed_plot.obs[key].astype(str).nunique() < 2:
            print(f"[skip] panels {key}: missing or constant", flush=True)
            continue
        groups = embed_plot.obs[key].astype("category")
        labels = [c for c in groups.cat.categories if (groups == c).any()]
        ncols = 4
        nrows = int(np.ceil(len(labels) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.4))
        axes = np.atleast_1d(axes).flatten()
        for ax, label in zip(axes, labels):
            n = int((groups == label).sum())
            sc.pl.umap(embed_plot, color=key, groups=[label], title=f"{label} ({n:,})",
                       ax=ax, show=False, legend_loc="none", na_in_legend=False, size=4)
            ax.set_xlabel("")
            ax.set_ylabel("")
        for ax in axes[len(labels):]:
            ax.set_visible(False)
        plt.tight_layout()
        savefig(fig_dir, run_id, f"umap_per_group_{key}", fig=fig)
    del embed_plot

    # Per-dimension statistics, with and without the vanished ones.
    drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, show=False)
    savefig(fig_dir, run_id, "latent_dimension_stats")
    drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, remove_vanished=True, show=False)
    savefig(fig_dir, run_id, "latent_dimension_stats_rmVanished")

    # Every non-vanished dimension on the UMAP (remove_vanished defaults to True).
    drvi.utils.pl.plot_latent_dims_in_umap(embed, show=False)
    savefig(fig_dir, run_id, "latent_dims_in_umap")

    # Which dimension responds to which group, sorted so a group's own dimensions sit together.
    for key in heatmap_keys:
        drvi.utils.pl.plot_latent_dims_in_heatmap(embed, key, title_col="title",
                                                  sort_by_categorical=True, show=False)
        savefig(fig_dir, run_id, f"latent_dims_in_heatmap_{key}_sorted")

    # Interpretability: OOD from decoder reconstructions, IND averaged over all cells.
    model.plot_interpretability_scores(embed, adata, show=False)
    savefig(fig_dir, run_id, "ood_interpretability_scores")
    model.plot_interpretability_scores(embed, adata, key="OOD_max_possible", show=False)
    savefig(fig_dir, run_id, "ood_max_interpretability_scores")
    model.plot_interpretability_scores(embed, adata, key="OOD_min_possible", show=False)
    savefig(fig_dir, run_id, "ood_min_interpretability_scores")
    model.plot_interpretability_scores(embed, adata, key="IND_linear_weighted_mean", show=False)
    savefig(fig_dir, run_id, "ind_linear_weighted_mean")


def main():
    args = parse_args()

    if not args.data_dir:
        sys.exit("set DATA_DIR (or pass --data-dir) to the directory holding the datasets")
    os.environ["DATA_DIR"] = str(args.data_dir)   # read by cell_set

    C.banner("03_3 DRVI run")

    run_id = C.run_id()  # drvi_tum_64_nomt

    phase_dir = Path(__file__).resolve().parent.parent   # 03_drvi_tumoral_epi/
    tum_dir = C.tum_dir()                                # $DATA_DIR/03_tum

    input_h5ad = C.hvg_path(".h5ad")
    # All genes: gets the latent space and the Leiden columns at the end.
    full_h5ad = C.path(".h5ad")

    model_path = tum_dir / f"model_{run_id}.pt"          # one file (see save_model)
    latent_h5ad = tum_dir / f"embed_{run_id}.h5ad"       # latent space + interpretability
    downstream_h5ad = C.path(f"_{run_id}.h5ad")          # 03_2 object + obsm['X_drvi']

    fig_dir = Path(args.fig_dir) if args.fig_dir else phase_dir / "figures" / f"03_3_{run_id}"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"run id      {run_id}", flush=True)
    print(f"input       {input_h5ad}")
    print(f"model       {model_path}")
    print(f"latent      {latent_h5ad}")
    print(f"downstream  {downstream_h5ad}")
    print(f"figures     {fig_dir}")
    print(f"scvi-tools {scvi.__version__}, DRVI {drvi.__version__}", flush=True)

    if not input_h5ad.exists():
        sys.exit(f"missing {input_h5ad}: run 03_2_subsetting/subsetting_all.sh first")

    # Checked before training rather than hours later: the downstream object needs it.
    if not full_h5ad.exists():
        sys.exit(f"missing {full_h5ad}: run 03_2_subsetting/subsetting_all.sh first")

    sc.settings.set_figure_params(dpi=300, facecolor="white")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(input_h5ad)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} HVGs", flush=True)

    expected = C.compartment()
    assert adata.obs[C.COMPARTMENT_KEY].astype(str).nunique() == 1, \
        "not a single-compartment object"
    assert adata.obs[C.COMPARTMENT_KEY].astype(str).iloc[0] == expected, \
        f"compartment is not {expected!r}: this is not the malignant object"
    assert "counts" in adata.layers, "raw counts layer missing, DRVI has nothing to train on"
    assert (adata.obs[C.STATUS_KEY].astype(str) != "not_tested").all(), \
        "not_tested cells in the input: 03_2 should have dropped them"
    print("compartment :", adata.obs[C.COMPARTMENT_KEY].astype(str).unique().tolist())
    print("cnv_status  :", sorted(adata.obs[C.STATUS_KEY].astype(str).unique().tolist()))
    print("cell_type   :", sorted(adata.obs[LABEL_KEY].astype(str).unique().tolist()))
    print(f"{adata.obs[BATCH_KEY].nunique()} cohorts", flush=True)

    DRVI.setup_anndata(adata, layer="counts", batch_key=BATCH_KEY)
    scvi.settings.seed = args.seed

    model = DRVI(
        adata,
        n_latent=C.N_LATENT,
        encoder_dims=[128, 128],
        decoder_dims=[128, 128],
        dispersion="gene-batch",
    )

    if args.overwrite or not model_path.exists():
        print(f">>> training {run_id}: {args.epochs} epochs, early stopping "
              f"{'on (patience 50)' if args.early_stopping else 'off'}", flush=True)
        model.train(
            max_epochs=args.epochs,
            early_stopping=args.early_stopping,
            early_stopping_patience=50,
        )
        save_model(model, model_path)
        print(f"[model] saved to {model_path}", flush=True)
    else:
        print(f"[have] {model_path} already exists, not retraining (--overwrite to force)",
              flush=True)

    # Reloaded, so a resumed run starts from the same state as a fresh one.
    model = load_model(model_path, adata)

    # Latent space as .X, one var per dimension with its stats and scores; read by 03_4.
    if args.overwrite or not latent_h5ad.exists():
        embed = ad.AnnData(model.get_latent_representation(), obs=adata.obs)

        print("Setting latent dimension stats ...", flush=True)
        model.set_latent_dimension_stats(embed, vanished_threshold=0.5)

        print("Calculating gene scores per factor ...", flush=True)
        model.calculate_interpretability_scores(embed, "OOD")
        model.calculate_interpretability_scores(embed, "IND")

        # Vanished dimensions are out from here on: they carry no signal, only noise.
        live_umap(embed)

        print(f"[write] {latent_h5ad}", flush=True)
        embed.write_h5ad(latent_h5ad)
    else:
        print(f"[have] {latent_h5ad}, reading it back", flush=True)
        embed = sc.read_h5ad(latent_h5ad)

    # No vanished dimension would mean 64 is too tight.
    n_vanished = int(embed.var["vanished"].sum())
    print(f"{n_vanished} vanished / {embed.n_vars} latent dimensions "
          f"({embed.n_vars - n_vanished} effectively used)", flush=True)
    if n_vanished == 0:
        print(f"[note] nothing vanished at n_latent={C.N_LATENT}: the space may be too tight",
              flush=True)

    # Stored in the embedding, so 03_4 - 03_9 read the same partition.
    leiden_keys = [C.leiden_key(r) for r in C.LEIDEN_RESOLUTIONS] + [LEIDEN_KEY]
    table_dir = phase_dir / "tables" / f"03_3_{run_id}"
    table_dir.mkdir(parents=True, exist_ok=True)
    sizes_csv = table_dir / f"drvi_leiden_clusters_{run_id}.csv"
    # Recomputed when a Leiden setting changed since the columns were written (embed.uns).
    if args.overwrite or not leiden_is_current(embed):
        sizes = drvi_leiden(embed)
        sizes.to_csv(sizes_csv, index=False)
        print(f"[table] {sizes_csv}", flush=True)
        print(f"[write] {latent_h5ad} (+ Leiden columns)", flush=True)
        embed.write_h5ad(latent_h5ad)
    else:
        print(f"[have] Leiden columns in {latent_h5ad.name}", flush=True)
    heatmap_keys = groupable(embed, HEATMAP_KEYS)

    # The 03_2 object (all genes) + obsm['X_drvi'] + Leiden, for steps needing genes and latent space.
    if args.overwrite or not downstream_h5ad.exists():
        print(f"[read] {full_h5ad}", flush=True)
        full = sc.read_h5ad(full_h5ad)

        # Realigned by cell name, not by order.
        assert set(embed.obs_names) == set(full.obs_names), \
            "the embedding and the 03_2 object disagree on cells"
        full.obsm["X_drvi"] = np.asarray(embed[full.obs_names].X, dtype=np.float32)
        assert np.isfinite(full.obsm["X_drvi"]).all()
        for key in leiden_keys:
            full.obs[key] = embed.obs.loc[full.obs_names, key].values

        print(f"Writing {downstream_h5ad.name} (gzip, a few minutes) ...", flush=True)
        full.write_h5ad(downstream_h5ad, compression="gzip")
        print(f"[write] {downstream_h5ad} "
              f"({full.n_obs:,} x {full.n_vars:,}, obsm['X_drvi'] "
              f"{full.obsm['X_drvi'].shape}, "
              f"{downstream_h5ad.stat().st_size / 1024 ** 3:.2f} GB on disk)", flush=True)
        del full
    else:
        print(f"[have] {downstream_h5ad} (--overwrite to rewrite)", flush=True)

    print(">>> figures", flush=True)
    try:
        leiden_resolution_grid(embed, fig_dir, run_id)
        draw_figures(model, adata, embed, fig_dir, run_id, heatmap_keys)
    except Exception:
        # Every output is already written: report and let the notebook redraw.
        print("[warn] the figures section failed; the artifacts above are complete",
              flush=True)
        traceback.print_exc()

    print(f"[ok] {run_id} -> {latent_h5ad}", flush=True)
    print("next: 03_4_signatures on that embedding", flush=True)


if __name__ == "__main__":
    main()
