#!/usr/bin/env python3
"""Library for route_a_null_tum.py: the size- and expression-matched random gene-set null.

Random sets are drawn bin for bin against the all-genes object, scored with
`sc.tl.score_genes` and z-scored within cohort as 03_6 does, then cached in $DATA_DIR/03_tum.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import cell_set as CS                    # noqa: E402
import signature_common as C             # noqa: E402

N_BINS = 25              # sc.tl.score_genes' own default, as 03_6 states it


def z_within(df: pd.DataFrame, cohorts: pd.Series) -> pd.DataFrame:
    """Z-score each column within `cohort`, 03_6's GROUPBY and 03_6's zero-variance rule."""
    g = df.groupby(cohorts.values, observed=True)
    z = g.transform(lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else s * 0.0)
    return z.fillna(0.0)


def expression_bins(adata, n_bins: int = N_BINS) -> pd.Series:
    """Every gene's expression bin, as sc.tl.score_genes builds its control bins (not exposed by scanpy)."""
    X = adata.X
    mean = np.asarray(X.mean(axis=0)).ravel()
    ranked = pd.Series(rankdata(mean, method="min"), index=adata.var_names)
    return pd.Series(pd.cut(ranked, n_bins, labels=False), index=adata.var_names)


def draw_matched_sets(genes: list[str], bins: pd.Series, n_sets: int, rng) -> list[list[str]]:
    """`n_sets` gene sets of the same size as `genes`, one gene per bin of each real gene,
    so expression level is matched as well as size."""
    by_bin = {b: np.asarray(idx) for b, idx in bins.groupby(bins).groups.items()}
    present = [g for g in genes if g in bins.index]
    wanted = bins.loc[present].values
    banned = set(present)
    sets = []
    for _ in range(n_sets):
        pick = []
        for b in wanted:
            pool = by_bin.get(b)
            if pool is None or len(pool) == 0:
                continue
            for _try in range(20):
                g = str(rng.choice(pool))
                if g not in banned and g not in pick:
                    pick.append(g)
                    break
        if len(pick) >= C.MIN_SIGNATURE_GENES:
            sets.append(pick)
    return sets


def build_random_scores(coll, sets: dict[str, list[str]], obs_names, cohorts,
                        n_per_sig: int, cache, overwrite: bool, seed: int) -> pd.DataFrame:
    """Per-cell z scores of the matched random sets, columns <SIGNATURE>__rnd<i>; cached (slow)."""
    if cache.exists() and not overwrite:
        print(f"[have] {cache.name}, reading it back", flush=True)
        df = pd.read_csv(cache, index_col=0)
        df.index = df.index.astype(str)
        return df.reindex(obs_names)

    full = CS.path(".h5ad")
    if not full.exists():
        print(f"[skip] {full.name} is not on disk, so the matched random null cannot be "
              f"built; the permutation null still is", flush=True)
        return pd.DataFrame(index=obs_names)

    print(f"[read] {full}  (the all-genes object, for the null only)", flush=True)
    adata = sc.read_h5ad(full)
    adata = adata[obs_names].copy()
    bins = expression_bins(adata)
    rng = np.random.default_rng(seed)

    raw = {}
    for name, genes in sets.items():
        if len(genes) < C.MIN_SIGNATURE_GENES:
            continue
        for i, pick in enumerate(draw_matched_sets(genes, bins, n_per_sig, rng)):
            key = f"{name}__rnd{i}"
            sc.tl.score_genes(adata, gene_list=pick, ctrl_size=len(pick), n_bins=N_BINS,
                              random_state=seed, use_raw=False, score_name=key)
            raw[key] = adata.obs[key].astype(float).values
        print(f"  {name:28s} {len(genes):4d} genes -> {n_per_sig} matched random sets",
              flush=True)

    if not raw:
        return pd.DataFrame(index=obs_names)
    df = z_within(pd.DataFrame(raw, index=adata.obs_names.astype(str)), cohorts)
    df.to_csv(cache)
    print(f"[write] {cache}  ({df.shape[0]} x {df.shape[1]})", flush=True)
    return df.reindex(obs_names)
