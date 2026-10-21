"""
02_2 integration: method implementations for the Python side.

Every method has the signature
    method(adata, batch_key, label_key=None, seed=0, model_dir=None,
           model_prefix="", retrain=False) -> AnnData
and returns what the metrics need for its output type:
    knn   (bbknn)                 -> corrected neighbour graph in .obsp / .uns
    embed (harmony, scvi, scanvi) -> obsm['X_emb']
    full + embed (scanorama)      -> corrected .X and obsm['X_emb']

BBKNN and Harmony wrap scib.integration. Scanorama is reimplemented because scib's
wrapper needs AnnData.concatenate (removed in anndata 0.13). scVI and scANVI copy scib
exactly but add a checkpoint in model_dir, so a crash after training is cheap.

The caller must import scib_compat before scib.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import scib

from model_paths import has_saved_model, saved_model


def bbknn(adata, batch_key, label_key=None, seed=0, model_dir=None,
          model_prefix="", retrain=False):
    """BBKNN: batch-balanced kNN graph. Output type: knn.

    scib.integration.bbknn recomputes PCA and returns a copy with the corrected
    neighbour graph in .obsp / .uns['neighbors']. It switches to
    neighbors_within_batch=25 above 1e5 cells, which is our regime.
    """
    return scib.integration.bbknn(adata, batch_key)


def harmony(adata, batch_key, label_key=None, seed=0, model_dir=None,
            model_prefix="", retrain=False):
    """Harmony (harmony-pytorch). Output type: embed.

    scib.integration.harmony imports `from harmony import harmonize`.
    It computes its own PCA on .X and writes the corrected
    embedding to obsm['X_emb']; .X is left untouched. It ignores any precomputed
    X_pca, so on the scaled variant the PCA it uses is the scaled one, as intended.
    """
    return scib.integration.harmony(adata, batch_key)


def _fit_or_load_scvi(adata, batch_key, model_dir, prefix, retrain, max_epochs=None):
    """Train an SCVI, or reload one from model_dir. Returns (model, net_adata).

    scib.integration.scvi's body plus a checkpoint; all result-affecting settings
    match scib (layers['counts'], n_latent=30, n_hidden=128, n_layers=2, nb,
    train_size=1.0, default epoch heuristic). setup_anndata runs on a copy.
    """
    from scvi.model import SCVI

    scib.utils.check_sanity(adata, batch_key, None)
    if "counts" not in adata.layers:
        raise TypeError("adata has no layers['counts']; scVI needs raw counts")

    net_adata = adata.copy()
    SCVI.setup_anndata(net_adata, layer="counts", batch_key=batch_key)

    if not retrain and has_saved_model(model_dir, prefix):
        print(f"[load] scVI model from {saved_model(model_dir, prefix)} "
              "(--retrain to ignore)", flush=True)
        return SCVI.load(model_dir, adata=net_adata, prefix=prefix), net_adata

    vae = SCVI(net_adata, gene_likelihood="nb", n_layers=2, n_latent=30, n_hidden=128)
    train_kwargs = {"train_size": 1.0}
    if max_epochs is not None:
        train_kwargs["max_epochs"] = max_epochs
    vae.train(**train_kwargs)

    if model_dir is not None:
        print(f"[save] scVI model -> {saved_model(model_dir, prefix)}", flush=True)
        vae.save(model_dir, prefix=prefix, overwrite=True)
    return vae, net_adata


def scvi(adata, batch_key, label_key=None, seed=0, model_dir=None,
         model_prefix="", retrain=False):
    """scVI on layers['counts'] (unscaled input). Output type: embed."""
    import scvi as scvi_tools

    scvi_tools.settings.seed = seed
    vae, _ = _fit_or_load_scvi(adata, batch_key, model_dir, model_prefix, retrain)
    adata.obsm["X_emb"] = vae.get_latent_representation()
    return adata


def scanvi(adata, batch_key, label_key=None, seed=0, model_dir=None,
           model_prefix="", retrain=False):
    """scANVI (needs label_key) on layers['counts']. Output type: embed.

    The scVI and scANVI stages are checkpointed separately, so a resume reloads
    whichever is on disk. Epochs follow scib.integration.scanvi (13 then 4 at 619k cells).
    """
    if label_key is None:
        raise ValueError("scanvi requires label_key (the cell-type annotation)")

    import scvi as scvi_tools
    from scvi.model import SCANVI

    scvi_tools.settings.seed = seed

    n_epochs_scvi = int(np.min([round((20000 / adata.n_obs) * 400), 400]))
    n_epochs_scanvi = int(np.min([10, np.max([2, round(n_epochs_scvi / 3.0)])]))

    # One prefix per stage, both in the same flat directory.
    scvi_prefix = f"{model_prefix}scvi_"
    scanvi_prefix = f"{model_prefix}scanvi_"

    # SCANVI.load rebuilds the full model, so the scVI stage is not needed.
    if not retrain and has_saved_model(model_dir, scanvi_prefix):
        print(f"[load] scANVI model from {saved_model(model_dir, scanvi_prefix)} "
              "(--retrain to ignore)", flush=True)
        net_adata = adata.copy()
        SCANVI.setup_anndata(net_adata, layer="counts", batch_key=batch_key,
                             labels_key=label_key, unlabeled_category="UnknownUnknown")
        scanvae = SCANVI.load(model_dir, adata=net_adata, prefix=scanvi_prefix)
        adata.obsm["X_emb"] = scanvae.get_latent_representation()
        return adata

    vae, _ = _fit_or_load_scvi(adata, batch_key, model_dir, scvi_prefix, retrain,
                               max_epochs=n_epochs_scvi)

    scanvae = SCANVI.from_scvi_model(
        scvi_model=vae,
        labels_key=label_key,
        unlabeled_category="UnknownUnknown",  # anything definitely not in the data
    )
    scanvae.train(max_epochs=n_epochs_scanvi, train_size=1.0)

    if model_dir is not None:
        print(f"[save] scANVI model -> {saved_model(model_dir, scanvi_prefix)}", flush=True)
        scanvae.save(model_dir, prefix=scanvi_prefix, overwrite=True)

    adata.obsm["X_emb"] = scanvae.get_latent_representation()
    return adata


def scanorama(adata, batch_key, label_key=None, seed=0, model_dir=None,
              model_prefix="", retrain=False):
    """Scanorama, reimplemented to preserve cell order. Output types: full + embed.

    Results are written back by original row position: corrected .X and obsm['X_emb'].
    """
    try:
        import scanorama
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "scanorama is not installed in this environment"
        ) from error

    if not isinstance(adata.obs[batch_key].dtype, pd.CategoricalDtype):
        adata.obs[batch_key] = adata.obs[batch_key].astype("category")
    batches = adata.obs[batch_key].cat.categories
    batch_values = adata.obs[batch_key].to_numpy()

    # Split by batch, remembering each batch's original row positions.
    splits = []
    row_index = []
    for batch in batches:
        rows = np.flatnonzero(batch_values == batch)
        if rows.size == 0:
            raise ValueError(f"batch {batch!r} is an empty category; drop it upstream")
        row_index.append(rows)
        splits.append(adata[rows].copy())

    # return_dimred=True also writes the joint embedding to obsm['X_scanorama'].
    corrected = scanorama.correct_scanpy(splits, return_dimred=True)

    emb_dim = corrected[0].obsm["X_scanorama"].shape[1]
    corrected_x = np.empty((adata.n_obs, adata.n_vars), dtype=np.float32)
    corrected_emb = np.empty((adata.n_obs, emb_dim), dtype=np.float32)
    written = np.zeros(adata.n_obs, dtype=bool)

    for rows, part in zip(row_index, corrected):
        block = part.X
        block = block.toarray() if sp.issparse(block) else np.asarray(block)
        corrected_x[rows] = block.astype(np.float32, copy=False)
        corrected_emb[rows] = part.obsm["X_scanorama"].astype(np.float32, copy=False)
        written[rows] = True

    assert written.all(), f"{int((~written).sum())} cells never corrected by scanorama"

    adata.X = corrected_x
    adata.obsm["X_emb"] = corrected_emb
    return adata


# Registry consumed by run_integration.py. label_key is only used by scanvi.
METHODS = {
    "bbknn": bbknn,
    "harmony": harmony,
    "scvi": scvi,
    "scanvi": scanvi,
    "scanorama": scanorama,
}
