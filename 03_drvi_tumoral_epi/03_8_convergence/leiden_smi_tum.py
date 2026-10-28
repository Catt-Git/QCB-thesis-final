#!/usr/bin/env python3
"""03_8, the DRVI-Leiden clusters: which DR direction owns which cluster (SMI).

SMI of every live direction against 03_3's Leiden clusters at --smi-resolution (1.5,
drvi_leiden), as in the DRVI tutorial; tables only (the figures are 03_3's). dr_naming_tum.py
lays the STATE_PICKS matches over the names, leiden_states_tum.py characterises the clusters.

Outputs: tables/03_8_convergence/smi_matrix_drvi_leiden_<res>_<run>.csv and
smi_matches_drvi_leiden_<res>_<run>.csv (pairs at or above --smi-threshold).

Usage:
    export DATA_DIR=<data root>
    python3 leiden_smi_tum.py
"""

from __future__ import annotations

import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402

CS = C.CS

STEP = "03_8_convergence"


def load_partitions(drvi_run: str) -> pd.DataFrame:
    """03_3's Leiden partitions on the live DRs: cells (by name) x resolution ('0.25', ...)."""
    import h5py
    from anndata.io import read_elem
    path = CS.tum_dir() / f"embed_{drvi_run}.h5ad"
    if not path.exists():
        sys.exit(f"missing {path}: run 03_3_drvi_run first")
    # Only obs: the partitions are written by 03_3, never recomputed here.
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
    keys = {str(r): CS.leiden_key(r) for r in CS.LEIDEN_RESOLUTIONS}
    missing = [k for k in keys.values() if k not in obs]
    if missing:
        sys.exit(f"{path.name} has no obs[{missing[0]!r}]: re-run 03_3_drvi_run/run_drvi_tum.py")
    parts = pd.DataFrame({r: obs[k].astype(int) for r, k in keys.items()})
    parts.index = obs.index.astype(str)
    return parts


def leiden_smi(hvg, args, drvi_run: str):
    """SMI of the live DR directions against 03_3's Leiden clusters at --smi-resolution."""
    tdir = C.table_dir(STEP)
    parts = load_partitions(drvi_run)
    embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad")
    embed.obs_names = embed.obs_names.astype(str)
    # An obs-only view on the DRVI input's cells carries the cluster labels to the benchmark.
    cells = hvg.obs_names.astype(str)
    view = ad.AnnData(obs=pd.DataFrame(index=cells))
    lab = parts[str(args.smi_resolution)].reindex(cells)
    view.obs[CS.leiden_key(args.smi_resolution)] = pd.Categorical(
        lab.astype(str), categories=[str(c) for c in sorted(lab.unique())])
    drvi_states_smi(embed, view, args.smi_resolution, args.smi_threshold, drvi_run, tdir)


def drvi_states_smi(embed, view, res: float, threshold: float, drvi_run: str, tdir):
    """SMI between every live DR direction and the clusters of one resolution, tables only.

    As the DRVI tutorial (directions clipped at zero, vanished ones dropped); pairs at or above
    `threshold` (0.7) are kept. The figures are the 03_3 notebook's.
    """
    from drvi.utils.metrics import DiscreteDisentanglementBenchmark

    key = f"drvi_leiden_{res}"
    if key not in view.obs:
        print(f"[skip] SMI: no {key} labels for --smi-resolution {res}", flush=True)
        return
    C.banner(f"SMI of the live DR directions against {key}")

    # Each side of a dimension is its own direction: + keeps x > 0, - keeps -x for x < 0.
    pos = embed[:, ~embed.var["vanished_positive_direction"].astype(bool)].copy()
    neg = embed[:, ~embed.var["vanished_negative_direction"].astype(bool)].copy()
    pos.var.index = pos.var["title"].astype(str) + "+"
    neg.var.index = neg.var["title"].astype(str) + "-"
    pos.X = np.asarray(pos.X).clip(min=0)
    neg.X = -np.asarray(neg.X).clip(max=0)
    directional = pd.concat([pos.to_df(), neg.to_df()], axis=1).loc[view.obs_names]
    print(f"    {directional.shape[1]} live directions x {view.obs[key].nunique()} clusters",
          flush=True)

    target = view.obs[key]
    bench = DiscreteDisentanglementBenchmark(directional.values,
                                             dim_titles=directional.columns,
                                             discrete_target=target, metrics=["SMI"],
                                             aggregation_methods=[])
    bench.evaluate()
    smi = bench.get_results_details()["SMI"]
    smi.index.name = "direction"
    # Clusters in partition order; the long table keeps only the pairs over the threshold.
    smi = smi[[c for c in target.cat.categories if c in smi.columns]]
    top = (smi.reset_index().melt(id_vars="direction", var_name="cluster", value_name="smi")
           .query("smi >= @threshold").sort_values("smi", ascending=False)
           .reset_index(drop=True))

    smi.to_csv(tdir / f"smi_matrix_{key}_{drvi_run}.csv")
    top.to_csv(tdir / f"smi_matches_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'smi_matrix_{key}_{drvi_run}.csv'}")
    print(f"[table] {tdir / f'smi_matches_{key}_{drvi_run}.csv'}  ({len(top)} pairs >= "
          f"{threshold})")
    best = smi.max(axis=0)
    print("    best SMI per cluster: " + ", ".join(f"{c} {v:.2f}" for c, v in best.items()),
          flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smi-resolution", type=float, default=CS.LEIDEN_RESOLUTION,
                   help="the 03_3 Leiden resolution whose clusters get the SMI screening "
                        "against the DR directions (default %(default)s); one of "
                        f"{list(CS.LEIDEN_RESOLUTIONS)}")
    p.add_argument("--smi-threshold", type=float, default=0.7,
                   help="SMI at or above which a direction-cluster pair is kept (default "
                        "%(default)s, as SMI_THRESHOLD in the 03_3 notebook)")
    return p.parse_args()


def main():
    args = parse_args()
    drvi_run = CS.run_id()
    CS.banner(f"{STEP} - SMI of the DR directions against the DRVI-Leiden clusters ({drvi_run})")
    # Backed: only the cell names are needed.
    hvg = ad.read_h5ad(C.HVG_H5AD, backed="r")
    leiden_smi(hvg, args, drvi_run)
    print("\ndone.")


if __name__ == "__main__":
    main()
