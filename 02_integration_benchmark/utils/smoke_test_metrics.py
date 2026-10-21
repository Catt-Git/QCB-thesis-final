"""
02 utils - smoke test for the metrics environment

Checks on a tiny sample, before the full grid, that scib, rpy2/anndata2ri and
all 12 metrics work for 'full' and 'embed'. The integration is the identity, so
values are meaningless: each metric only has to return a number.

Input : $DATA_DIR/shiao.h5ad                     (first run only)
        $DATA_DIR/shiao_hvg_2k_unintegrated_list.csv  (first run only)
Output: $SMOKE_DIR/smoke_fixture.h5ad            (sample, reused afterwards)
        $SMOKE_DIR/smoke_test_metrics.csv        (values obtained)
        $SMOKE_DIR/smoke_test_{full,embed}_nmi.txt  (scib NMI-vs-resolution scan)
SMOKE_DIR defaults to utils/smoke_test/ (disposable, not in git).

Build the fixture once locally, copy it to the cluster, run the test there.

Usage:
  export DATA_DIR=<data root>
  python3 smoke_test_metrics.py    # locally, once: builds the fixture
  scp smoke_test/smoke_fixture.h5ad <cluster>:<repo>/02_integration_benchmark/utils/smoke_test/
  sbatch submit_smoke_test.slurm   # on the cluster, from 02_integration_benchmark/utils

Exit code 0 if every expected metric produced a finite value, 1 otherwise.
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scib_compat  # noqa: F401  (restores the removed APIs: must precede scib)

import scanpy as sc
import scib


# CONFIG

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "shiao.h5ad")
HVG_PATH = os.path.join(DATA_DIR, "shiao_hvg_2k_unintegrated_list.csv")
SMOKE_DIR = os.environ.get(
    "SMOKE_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "smoke_test")
)
os.makedirs(SMOKE_DIR, exist_ok=True)
FIXTURE_PATH = os.path.join(SMOKE_DIR, "smoke_fixture.h5ad")

BATCH_KEY = "cohort"
LABEL_KEY = "cell_type"
ORGANISM = "human"
SEED = 0

# Sample composition; rare labels are needed for the isolated label metrics.
N_COHORTS = 4
N_TOP_LABELS = 5
N_RARE_LABELS = 2
MIN_CELLS_RARE_LABEL = 100
MAX_CELLS_PER_GROUP = 250

# Expected metrics per output type (HVG conservation only for 'full').
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
}

OUT_CSV = os.path.join(SMOKE_DIR, "smoke_test_metrics.csv")


def report_environment():
    """Print the stack versions and check that R is reachable through rpy2."""
    import importlib.metadata as md

    print("=" * 70, flush=True)
    print("ENVIRONMENT", flush=True)
    print("=" * 70, flush=True)
    print(f"  python      {sys.version.split()[0]}  ({sys.executable})", flush=True)
    for pkg in ["scanpy", "anndata", "scib", "rpy2", "anndata2ri", "numpy", "pandas"]:
        try:
            print(f"  {pkg:<11} {md.version(pkg)}", flush=True)
        except md.PackageNotFoundError:
            print(f"  {pkg:<11} NOT INSTALLED", flush=True)

    # scib_compat starts an embedded R, so report its version.
    import rpy2.robjects as ro

    print(f"  R           {ro.r('R.version.string')[0]}", flush=True)
    print(flush=True)


def build_fixture():
    """Small sample stratified on (cohort, cell_type), with a few rare cell types."""
    print(f"[fixture] reading {IN_PATH} (once only)...", flush=True)
    adata = sc.read_h5ad(IN_PATH)
    print(f"[fixture] full object: {adata.n_obs} x {adata.n_vars}", flush=True)

    cohorts = adata.obs[BATCH_KEY].value_counts().index[:N_COHORTS].tolist()
    adata = adata[adata.obs[BATCH_KEY].isin(cohorts)].copy()

    counts = adata.obs[LABEL_KEY].value_counts()
    top_labels = counts.index[:N_TOP_LABELS].tolist()
    rare_pool = counts[counts >= MIN_CELLS_RARE_LABEL]
    rare_labels = [lb for lb in rare_pool.index[::-1] if lb not in top_labels][:N_RARE_LABELS]
    labels = top_labels + rare_labels
    adata = adata[adata.obs[LABEL_KEY].isin(labels)].copy()

    print(f"[fixture] selected cohorts : {cohorts}", flush=True)
    print(f"[fixture] abundant labels  : {top_labels}", flush=True)
    print(f"[fixture] rare labels      : {rare_labels}", flush=True)

    # Stratified subsampling, with a fixed seed for reproducibility.
    rng = np.random.default_rng(SEED)
    keep = []
    for _, idx in adata.obs.groupby([BATCH_KEY, LABEL_KEY], observed=True).groups.items():
        idx = np.asarray(idx)
        if len(idx) > MAX_CELLS_PER_GROUP:
            idx = rng.choice(idx, MAX_CELLS_PER_GROUP, replace=False)
        keep.append(idx)
    adata = adata[np.concatenate(keep)].copy()

    # Drop empty categories (scib would return spurious NaNs).
    for key in [BATCH_KEY, LABEL_KEY]:
        adata.obs[key] = adata.obs[key].cat.remove_unused_categories()

    # Genes zeroed out by the subsampling break PCA and HVG overlap.
    sc.pp.filter_genes(adata, min_cells=1)

    # Graph and PCA of the full object are invalid on the subset.
    for slot in ["neighbors", "pca", "umap"]:
        adata.uns.pop(slot, None)
    for slot in ["distances", "connectivities"]:
        adata.obsp.pop(slot, None)

    print(f"[fixture] sample: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{adata.obs[BATCH_KEY].nunique()} batches, "
          f"{adata.obs[LABEL_KEY].nunique()} labels", flush=True)

    adata.write_h5ad(FIXTURE_PATH, compression="gzip")
    print(f"[fixture] wrote {FIXTURE_PATH}", flush=True)
    return adata


def load_hvg(adata):
    """HVGs from the csv if available, otherwise from the flag already in .var."""
    if os.path.exists(HVG_PATH):
        hvg = pd.read_csv(HVG_PATH, header=None)[0].astype(str).tolist()
        hvg = [g for g in hvg if g in adata.var_names]
        if hvg:
            return hvg
    return adata.var_names[adata.var["highly_variable"]].tolist()


def run_metrics(adata_pre, hvg, type_):
    """Reproduce the same sequence as 02_4_metrics on an identity integration."""
    print("=" * 70, flush=True)
    print(f"METRICS  --type {type_}", flush=True)
    print("=" * 70, flush=True)

    embed = "X_pca" if type_ == "full" else "X_emb"

    adata_int = adata_pre[:, hvg].copy()
    sc.pp.pca(adata_int, n_comps=30, svd_solver="arpack", random_state=SEED)
    if type_ == "embed":
        adata_int.obsm["X_emb"] = adata_int.obsm["X_pca"].copy()

    # Same reduce_data options as 02_4.
    scib.preprocessing.reduce_data(
        adata_int, n_top_genes=None, neighbors=True, use_rep=embed, pca=True, umap=False
    )

    results = scib.me.metrics(
        adata_pre,
        adata_int,
        verbose=True,
        batch_key=BATCH_KEY,
        label_key=LABEL_KEY,
        embed=embed,
        type_=type_,
        cluster_nmi=os.path.join(SMOKE_DIR, f"smoke_test_{type_}_nmi.txt"),
        hvg_score_=(type_ == "full"),
        silhouette_=True,
        nmi_=True,
        nmi_method="arithmetic",
        nmi_dir=None,
        ari_=True,
        pcr_=True,
        cell_cycle_=True,
        organism=ORGANISM,
        isolated_labels_=True,
        n_isolated=None,
        graph_conn_=True,
        lisi_graph_=True,
        trajectory_=False,
    )
    return results.iloc[:, 0]


def main():
    report_environment()

    if os.path.exists(FIXTURE_PATH):
        print(f"[fixture] reusing {FIXTURE_PATH}", flush=True)
        adata_pre = sc.read_h5ad(FIXTURE_PATH)
        print(f"[fixture] sample: {adata_pre.n_obs} x {adata_pre.n_vars}", flush=True)
    else:
        adata_pre = build_fixture()

    hvg = load_hvg(adata_pre)
    print(f"[fixture] HVGs available in the sample: {len(hvg)}\n", flush=True)

    # The reference needs its own PCA and graph.
    sc.pp.pca(adata_pre, n_comps=30, svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(adata_pre, random_state=SEED)

    obtained = {}
    for type_ in ["full", "embed"]:
        obtained[type_] = run_metrics(adata_pre.copy(), hvg, type_)


    # OUTCOME
    print("\n" + "=" * 70, flush=True)
    print("OUTCOME", flush=True)
    print("=" * 70, flush=True)

    failed = []
    rows = {}
    for type_, expected in EXPECTED.items():
        values = obtained[type_]
        rows[type_] = values
        print(f"\n--type {type_}", flush=True)
        for metric in expected:
            value = values.get(metric, np.nan)
            ok = np.isfinite(value) if isinstance(value, (int, float, np.floating)) else False
            print(f"  {'ok ' if ok else 'NaN'}  {metric:<28} {value}", flush=True)
            if not ok:
                failed.append(f"{type_}/{metric}")

    pd.DataFrame(rows).to_csv(OUT_CSV)
    print(f"\nValues saved to {OUT_CSV}", flush=True)

    if failed:
        print(f"\nSMOKE TEST FAILED: {len(failed)} metrics without a value", flush=True)
        for name in failed:
            print(f"  - {name}", flush=True)
        print(
            "\nA NaN metric here is not a biological result: it means the "
            "environment is unable to compute it. Fix this before launching the "
            "full grid.",
            flush=True,
        )
        return 1

    print("\nSMOKE TEST PASSED: every expected metric produced a value.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
