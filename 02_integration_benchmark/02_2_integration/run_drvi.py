#!/usr/bin/env python
"""
02_2 integration: DRVI, headless (the compute half of shiao_drvi_128.ipynb).

The latent size is chosen by eye in the notebook; this trains a chosen size
without a kernel (terminal or submit_drvi.slurm), with the notebook's exact
settings. Outputs under DATA_DIR, run id drvi_unscaled_<n_latent>:

    02_drvi/model_drvi_<N>.pt        the trained model
    02_drvi/embed_<run_id>.h5ad      latent space + dimension stats + OOD/IND scores
    02_integration/<run_id>.h5ad     benchmark output, obsm['X_emb'] (scored by 02_4)
    02_embeddings/<run_id>.npy       the embedding alone, for 02_3

The notebook's figures go to figures/<run_id>/ (a plotting failure does not fail
the run). Existing outputs are reused; --overwrite recomputes everything.

Usage:
  export DATA_DIR=<data root>
  python3 run_drvi.py                    # n_latent 128, the run in the grid
  python3 run_drvi.py --n-latent 64      # another size, side by side with the first
  python3 run_drvi.py --overwrite        # retrain and rewrite everything
  python3 run_drvi.py --early-stopping   # stop early instead (see --epochs)
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: figures are written to disk, never shown
import matplotlib.pyplot as plt  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402
import scvi  # noqa: E402
import drvi  # noqa: E402
from drvi.model import DRVI  # noqa: E402

BATCH_KEY = "cohort"       # the batch to correct, as everywhere else in the benchmark
LABEL_KEY = "cell_type"

# Same UMAP panels and names as the notebook.
UMAP_QC_KEYS = {
    "cell_type": "cell_type",
    "fraction": "fraction",
    "cohort": "cohort",
    "treatment": "treatment",
    "response": "response",
    "phase": "phase",
    "n_genes_by_counts": "n_genes_by_counts",
    "total_counts": "total_counts",
    "mito": "pct_counts_mt",
    "ribo": "pct_counts_ribo",
    "size_factors": "size_factors",
}
UMAP_COMBINED_KEYS = ["cell_type", "fraction", "cohort", "treatment", "response", "phase"]
UMAP_SEED = 0


def parse_args():
    p = argparse.ArgumentParser(description="Run DRVI on the prepared Shiao object")
    p.add_argument("-n", "--n-latent", type=int, default=128,
                   help="latent dimensions; the run id follows it [default: 128]")
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--epochs", type=int, default=400,
                   help="epochs to train; scvi-tools warms the KL term up over exactly "
                        "this many by default [default: 400]")
    p.add_argument("--early-stopping", action="store_true",
                   help="stop after 50 epochs without improvement. Off by default: the KL "
                        "warmup lasts --epochs epochs, so a run that stops earlier never "
                        "trains at full KL weight, and it is that KL pressure that makes "
                        "the unused latent dimensions vanish. Pair it with a shorter warmup "
                        "if you turn it on")
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR"),
                   help="directory holding the prepared inputs [default: $DATA_DIR]")
    p.add_argument("--fig-dir", default=None,
                   help="where the figures go [default: the repo's "
                        "02_integration_benchmark/figures/<run_id>]")
    p.add_argument("--overwrite", action="store_true",
                   help="retrain and rewrite everything instead of reusing what is on disk")
    return p.parse_args()


def save_model(model, path):
    """Save a DRVI model as the single file `path` (save() only writes <dir>/model.pt)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent) as scratch:
        model.save(scratch, overwrite=True)
        Path(scratch, "model.pt").replace(path)
    return path


def load_model(path, adata):
    """Load a model saved by save_model, via a scratch dir symlinking model.pt."""
    path = Path(path).resolve()
    with tempfile.TemporaryDirectory(dir=path.parent) as scratch:
        os.symlink(path, Path(scratch, "model.pt"))
        return DRVI.load(scratch, adata)


def savefig(fig_dir, run_id, name, fig=None, dpi=300):
    """Save and close a pyplot figure (current by default) as <name>_<run_id>.png."""
    fig = plt.gcf() if fig is None else fig
    path = Path(fig_dir) / f"{name}_{run_id}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path}", flush=True)
    return path


def draw_figures(model, adata, embed, fig_dir, run_id):
    """The notebook's plotting section. Cross-method QC panels are left to 02_3."""
    sc.settings.figdir = fig_dir

    # Same seeded point order as the notebook.
    order = np.random.default_rng(UMAP_SEED).permutation(embed.n_obs)
    embed_plot = ad.AnnData(
        obs=embed.obs.iloc[order].copy(),
        obsm={"X_umap": embed.obsm["X_umap"][order]},
        # Palettes carried over from the input object, so a cell type keeps the same
        # colour.
        uns={k: v for k, v in adata.uns.items() if k.endswith("_colors")},
    )

    for label, col in UMAP_QC_KEYS.items():
        if col not in embed_plot.obs:
            print(f"[skip] umap {label}: no obs column {col!r}", flush=True)
            continue
        sc.pl.umap(embed_plot, color=col, show=False, save=f"_{label}_{run_id}.png")
        plt.close("all")

    combined = [k for k in UMAP_COMBINED_KEYS if k in embed_plot.obs]
    with plt.rc_context({"figure.figsize": (7, 7)}):
        sc.pl.umap(embed_plot, color=combined, ncols=2, wspace=0.8,
                   hspace=0.25, show=False, save=f"_combined_{run_id}.png")
        plt.close("all")

    # One panel per cell type: ~50 labels are unreadable in a single UMAP.
    cell_types = [c for c in embed_plot.obs["cell_type"].cat.categories
                  if (embed_plot.obs["cell_type"] == c).any()]
    ncols = 4
    nrows = int(np.ceil(len(cell_types) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.4))
    axes = np.atleast_1d(axes).flatten()
    for ax, ct in zip(axes, cell_types):
        sc.pl.umap(embed_plot, color="cell_type", groups=[ct], title=ct, ax=ax,
                   show=False, legend_loc="none", na_in_legend=False, size=4)
        ax.set_xlabel("")
        ax.set_ylabel("")
    for ax in axes[len(cell_types):]:
        ax.set_visible(False)
    plt.tight_layout()
    savefig(fig_dir, run_id, "umap_per_cell_type", fig=fig)
    del embed_plot

    # Per-dimension stats; the rmVanished one drives the choice of n_latent.
    drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, show=False)
    savefig(fig_dir, run_id, "latent_dimension_stats")
    drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, remove_vanished=True, show=False)
    savefig(fig_dir, run_id, "latent_dimension_stats_rmVanished")

    # Every non-vanished dimension on the UMAP (remove_vanished defaults to True).
    drvi.utils.pl.plot_latent_dims_in_umap(embed, show=False)
    savefig(fig_dir, run_id, "latent_dims_in_umap")

    # Dimension vs cell type, sorted by dimension and then by label.
    drvi.utils.pl.plot_latent_dims_in_heatmap(embed, "cell_type", title_col="title", show=False)
    savefig(fig_dir, run_id, "latent_dims_in_heatmap_cell_type")
    drvi.utils.pl.plot_latent_dims_in_heatmap(embed, "cell_type", title_col="title",
                                              sort_by_categorical=True, show=False)
    savefig(fig_dir, run_id, "latent_dims_in_heatmap_cell_type_sorted")

    # The same heatmap against the batch, then treatment, then phase.
    for key in ("cohort", "treatment", "phase"):
        drvi.utils.pl.plot_latent_dims_in_heatmap(embed, key, title_col="title", show=False)
        savefig(fig_dir, run_id, f"latent_dims_in_heatmap_{key}")

    # Interpretability: OOD (decoder-based, specific genes), its min/max halves, IND (broader).
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
        sys.exit("set DATA_DIR (or pass --data-dir) to the directory holding the prepared inputs")

    # Latent size in the run id, so different sizes do not overwrite each other.
    run_id = f"drvi_unscaled_{args.n_latent}"

    script_dir = Path(__file__).resolve().parent
    phase_dir = script_dir.parent                     # 02_integration_benchmark/
    data_dir = Path(args.data_dir)

    # Same input as every other method.
    input_h5ad = data_dir / "shiao_hvg_2k.h5ad"

    # Outputs under DATA_DIR, not the repo (too large for GitHub).
    model_dir = data_dir / "02_drvi"                              # models and latent spaces, flat
    model_path = model_dir / f"model_drvi_{args.n_latent}.pt"     # one file (see save_model)
    latent_h5ad = model_dir / f"embed_{run_id}.h5ad"              # latent space + interpretability
    integrated_h5ad = data_dir / "02_integration" / f"{run_id}.h5ad"   # scored object, obsm['X_emb']
    emb_npy = data_dir / "02_embeddings" / f"{run_id}.npy"             # small durable copy

    fig_dir = Path(args.fig_dir) if args.fig_dir else phase_dir / "figures" / run_id

    for d in (model_dir, integrated_h5ad.parent, emb_npy.parent, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"run id      {run_id}", flush=True)
    print(f"input       {input_h5ad}")
    print(f"model       {model_path}")
    print(f"latent      {latent_h5ad}")
    print(f"integrated  {integrated_h5ad}")
    print(f"embedding   {emb_npy}")
    print(f"figures     {fig_dir}")
    print(f"scvi-tools {scvi.__version__}, DRVI {drvi.__version__}", flush=True)

    if not input_h5ad.exists():
        sys.exit(f"missing {input_h5ad}: run 02_1_prepare first")

    sc.settings.set_figure_params(dpi=300, facecolor="white")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(input_h5ad)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} HVGs", flush=True)
    assert "counts" in adata.layers, "raw counts layer missing, DRVI has nothing to train on"

    # Setup the anndata for scvi-tools
    DRVI.setup_anndata(
        adata,
        layer="counts",       # raw unnormalised data (used by default)
        batch_key=BATCH_KEY,  # the batch to fix
    )

    # Setting seed
    scvi.settings.seed = args.seed

    # Setup the model
    model = DRVI(
        adata,
        n_latent=args.n_latent,
        encoder_dims=[256, 128],
        decoder_dims=[128, 256],
        dispersion="gene-batch",   # really important parameter for TNBC
    )

    # Train only if no saved model, or --overwrite (GPU if available).
    if args.overwrite or not model_path.exists():
        print(f">>> training {run_id}: {args.epochs} epochs, early stopping "
              f"{'on (patience 50)' if args.early_stopping else 'off'}", flush=True)
        model.train(
            max_epochs=args.epochs,
            early_stopping=args.early_stopping,  # off: the KL warmup needs all --epochs
            early_stopping_patience=50,          # only read when --early-stopping is given
        )
        save_model(model, model_path)
        print(f"[model] saved to {model_path}", flush=True)
    else:
        print(f"[have] {model_path} already exists, not retraining (--overwrite to force)",
              flush=True)

    # Reload, so fresh and resumed runs continue from the same state.
    model = load_model(model_path, adata)

    # Latent space as .X, one var per dimension with stats and scores.
    if args.overwrite or not latent_h5ad.exists():
        embed = ad.AnnData(model.get_latent_representation(), obs=adata.obs)

        # Latent dimension statistics
        print("Setting latent dimension stats ...", flush=True)
        model.set_latent_dimension_stats(embed, vanished_threshold=0.5)

        # Interpretability gene scores
        print("Calculating gene scores per factor ...", flush=True)
        # OOD: from decoder reconstructions (faster)
        model.calculate_interpretability_scores(embed, "OOD")
        # IND: iterates over all cells
        model.calculate_interpretability_scores(embed, "IND")

        print("Dimension reduction ...", flush=True)
        sc.pp.neighbors(embed, n_neighbors=15, use_rep="X", n_pcs=embed.X.shape[1])
        sc.tl.umap(embed)
        sc.pp.pca(embed)

        print(f"[write] {latent_h5ad}", flush=True)
        embed.write_h5ad(latent_h5ad)
    else:
        print(f"[have] {latent_h5ad}, reading it back", flush=True)
        embed = sc.read_h5ad(latent_h5ad)

    n_vanished = int(embed.var["vanished"].sum())
    print(f"{n_vanished} vanished / {embed.n_vars} latent dimensions "
          f"({embed.n_vars - n_vanished} effectively used)", flush=True)

    # Benchmark output, same shape as run_integration.py's embed methods.
    if args.overwrite or not integrated_h5ad.exists():
        # Realign on adata's cells: the metrics depend on cell order.
        assert set(embed.obs_names) == set(adata.obs_names), "embed and adata disagree on cells"

        integrated = adata.copy()
        integrated.obsm["X_emb"] = np.asarray(embed[adata.obs_names].X, dtype=np.float32)
        assert np.isfinite(integrated.obsm["X_emb"]).all()

        # The metrics never read raw counts; dropping them keeps the file light.
        if "counts" in integrated.layers:
            del integrated.layers["counts"]

        print(f"[write] {integrated_h5ad}", flush=True)
        integrated.write_h5ad(integrated_h5ad, compression="gzip")

        np.save(emb_npy, integrated.obsm["X_emb"])
        print(f"[write] {emb_npy}  {integrated.obsm['X_emb'].shape}", flush=True)

        del integrated
    else:
        print(f"[have] {integrated_h5ad} (--overwrite to rewrite)", flush=True)

    print(">>> figures", flush=True)
    try:
        draw_figures(model, adata, embed, fig_dir, run_id)
    except Exception:
        # Artifacts are already written; the notebook can redraw the broken plot.
        print("[warn] the figures section failed; the artifacts above are complete",
              flush=True)
        traceback.print_exc()

    print(f"[ok] {run_id} -> {integrated_h5ad}", flush=True)


if __name__ == "__main__":
    main()
