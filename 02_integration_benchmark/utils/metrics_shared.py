"""
02 utils - shared preparation used by 02_4_metrics' metrics.py

Everything metrics.py does before scib.me.metrics, so every job prepares data
the same way:
  1. load_and_align: same cells in the same order (hard error otherwise); batch
     and label columns copied from the reference.
  2. ensure_reference_reduced: PCA and graph on the reference only if missing.
  3. reduce_integrated: scib's reduce_data by output type, umap=False.

Imports scib_compat before scib, so callers need not.
"""

from __future__ import annotations

import functools
import inspect
import os
import sys

import numpy as np
import scanpy as sc

# scib_compat must be imported before scib.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scib_compat  # noqa: F401,E402  (must precede scib)
import scib  # noqa: E402

# scib's PCA parameters, for references without a stored PCA.
PCA_N_COMPS = 50
SEED = 0

# scib output types, and the metric flags each one enables.
RESULT_TYPES = ("full", "embed", "knn")

# The twelve metrics, in scib's order; keep in sync with merge_metrics.SCORED_METRICS.
ALL_METRICS = (
    "NMI_cluster/label", "ARI_cluster/label", "ASW_label", "ASW_label/batch",
    "PCR_batch", "cell_cycle_conservation", "isolated_label_F1",
    "isolated_label_silhouette", "graph_conn", "iLISI", "cLISI",
    "hvg_overlap",
)

# Metrics expected to be finite per output type; the rest are NaN by design.
EXPECTED = {
    "full": [
        "NMI_cluster/label", "ARI_cluster/label", "ASW_label", "ASW_label/batch",
        "PCR_batch", "cell_cycle_conservation", "isolated_label_F1",
        "isolated_label_silhouette", "graph_conn", "iLISI", "cLISI",
        "hvg_overlap",
    ],
    "embed": [
        "NMI_cluster/label", "ARI_cluster/label", "ASW_label", "ASW_label/batch",
        "PCR_batch", "cell_cycle_conservation", "isolated_label_F1",
        "isolated_label_silhouette", "graph_conn", "iLISI", "cLISI",
    ],
    # knn (BBKNN): graph-based metrics only.
    "knn": [
        "NMI_cluster/label", "ARI_cluster/label", "isolated_label_F1",
        "graph_conn", "iLISI", "cLISI",
    ],
}


def use_igraph_leiden(n_iterations: int = 2, verbose: bool = False) -> bool:
    """Make scib's optimised clustering use igraph's leiden instead of leidenalg.

    leidenalg takes ~20 h here. Rebinding sc.tl.leiden is enough because scib
    looks it up at call time. Only NMI/ARI change slightly, so every run must use
    the same flavor.

    :returns: True if the swap happened, False on a scanpy without ``flavor``.
    """
    if "flavor" not in inspect.signature(sc.tl.leiden).parameters:
        print("[cluster] WARNING: this scanpy's leiden has no `flavor` argument; "
              "falling back to leidenalg. Expect the clustering to take hours.",
              flush=True)
        return False

    original = sc.tl.leiden
    if getattr(original, "_igraph_flavor", False):  # already patched
        return True

    @functools.wraps(original)  # keeps __name__ == "leiden", which scib prints
    def leiden_igraph(adata, **kwargs):
        kwargs.setdefault("flavor", "igraph")
        kwargs.setdefault("n_iterations", n_iterations)
        kwargs.setdefault("directed", False)  # required by the igraph flavor
        return original(adata, **kwargs)

    leiden_igraph._igraph_flavor = True
    sc.tl.leiden = leiden_igraph
    if verbose:
        print(f"[cluster] leiden flavor=igraph, n_iterations={n_iterations}", flush=True)
    return True


def metric_flags(type_: str, scaling: str = "unscaled") -> dict:
    """scib.me.metrics flags for one output type (scib-pipeline defaults, no trajectory).

    No PCR batch on scaled runs: the scaled reference has no batch means left, so the
    score is 0 by construction.
    """
    if type_ not in RESULT_TYPES:
        raise ValueError(f"unknown output type {type_!r}, expected one of {RESULT_TYPES}")

    flags = dict(
        silhouette_=True, nmi_=True, ari_=True, pcr_=True, cell_cycle_=True,
        isolated_labels_=True, hvg_score_=True, graph_conn_=True, lisi_graph_=True,
        trajectory_=False,
    )
    if type_ == "embed":
        flags["hvg_score_"] = False
    elif type_ == "knn":
        flags["silhouette_"] = False
        flags["pcr_"] = False
        flags["cell_cycle_"] = False
        flags["hvg_score_"] = False
    if scaling == "scaled":
        flags["pcr_"] = False
    return flags


def embed_key(type_: str) -> str:
    """The ``embed`` argument scib.me.metrics expects for one output type."""
    return "X_emb" if type_ == "embed" else "X_pca"


def load_and_align(uncorrected, integrated, batch_key, label_key, verbose=False):
    """Read reference and integrated object and enforce identical cell order.

    :returns: ``(adata, adata_int)`` with ``batch_key``/``label_key`` on
        ``adata_int`` guaranteed to match the reference cell for cell.
    :raises ValueError: if the two objects are not the same cells in the same
        order.
    """
    if verbose:
        print(f"[align] reference : {uncorrected}", flush=True)
    adata = sc.read_h5ad(uncorrected)
    if verbose:
        print(f"[align] integrated: {integrated}", flush=True)
    adata_int = sc.read_h5ad(integrated)

    if adata.n_obs != adata_int.n_obs:
        raise ValueError(
            f"cell count differs: reference has {adata.n_obs}, integrated has "
            f"{adata_int.n_obs}. The two objects must be the same cells."
        )

    ref_names = adata.obs_names.to_numpy()
    int_names = adata_int.obs_names.to_numpy()
    if not np.array_equal(ref_names, int_names):
        # Permuted vs different cells: both fatal, different upstream bugs.
        if set(ref_names) == set(int_names):
            first = int(np.argmax(ref_names != int_names))
            raise ValueError(
                "cell order differs between reference and integrated object: same "
                f"cells, permuted. First mismatch at position {first}: reference "
                f"{ref_names[first]!r} vs integrated {int_names[first]!r}. The "
                "metrics compare cells positionally; fix the integration output "
                "before scoring it (see check_integrations.py)."
            )
        missing = len(set(ref_names) - set(int_names))
        extra = len(set(int_names) - set(ref_names))
        raise ValueError(
            "reference and integrated object are not the same cells: "
            f"{missing} in reference only, {extra} in integrated only."
        )

    # Copy batch/label from the reference: some outputs drop or recategorise them.
    for key in (batch_key, label_key):
        if key not in adata.obs:
            raise ValueError(f"reference is missing obs[{key!r}]")
        adata_int.obs[key] = adata.obs[key].to_numpy()
        adata_int.obs[key] = adata_int.obs[key].astype("category")

    return adata, adata_int


def ensure_reference_reduced(adata, verbose=False):
    """Give the reference a PCA (with stored variance) and a graph, only if missing."""
    has_pca = "X_pca" in adata.obsm and "pca" in adata.uns
    if not has_pca:
        if verbose:
            print(f"[ref] recomputing PCA ({PCA_N_COMPS} comps, arpack)", flush=True)
        sc.pp.pca(adata, n_comps=PCA_N_COMPS, svd_solver="arpack", random_state=SEED)

    has_graph = "neighbors" in adata.uns and "connectivities" in adata.obsp
    if not has_graph:
        if verbose:
            print("[ref] recomputing neighbour graph", flush=True)
        sc.pp.neighbors(adata, random_state=SEED)

    return adata


def reduce_integrated(adata_int, type_, n_hvgs=None, verbose=False):
    """Run scib's reduce_data by type; knn keeps the method's own graph."""
    embed = embed_key(type_)
    if type_ == "knn":
        neighbors, pca = False, False
    else:
        neighbors, pca = True, True

    if verbose:
        print(f"[reduce] type={type_} use_rep={embed} neighbors={neighbors} "
              f"pca={pca} n_top_genes={n_hvgs}", flush=True)

    scib.preprocessing.reduce_data(
        adata_int,
        n_top_genes=n_hvgs,
        neighbors=neighbors,
        use_rep=embed,
        pca=pca,
        umap=False,
    )
    return adata_int
