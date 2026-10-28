#!/usr/bin/env python3
"""The null distribution of Route A's statistic, per collection: what `ROUTE_A_RHO_MIN` is against.

For one dimension-direction, the best rho over the collection's gene sets when each set is
replaced by a size- and expression-matched random one (scored and z-scored as in 03_6, cached
in $DATA_DIR/03_tum). Per collection, since the null moves with the number, overlap and size
of the sets.

Usage:
    export DATA_DIR=<data root>
    python3 route_a_null_tum.py                    # all three
    python3 route_a_null_tum.py --collection emt
    python3 route_a_null_tum.py --n-random 10      # a wider null
"""

from __future__ import annotations

import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "utils"))
import signature_presence_tum as SP  # noqa: E402
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402
import cell_set as CS  # noqa: E402

# Bars to price: 0.30 is ROUTE_A_RHO_MIN, the rest the round numbers around it.
BARS = (0.20, 0.225, 0.25, 0.275, 0.30, 0.325, 0.35, 0.372, 0.40)
PCTS = (50, 75, 90, 95, 97.5, 99)

# Pseudo-collections resampled from the cached draws, to smooth the maximum.
N_PSEUDO = 400


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--collection", choices=("scie", "emt", "gavish_tnbc"), action="append",
                   help="repeatable; default is all three")
    p.add_argument("--n-random", type=int, default=5,
                   help="matched random sets per gene set (default 5); a cache already on "
                        "disk is reused whatever this says, so raise it and delete the cache "
                        "to widen the null")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def ranked_z(X: np.ndarray) -> np.ndarray:
    """Columns ranked and standardised, so a dot product is Spearman."""
    R = np.apply_along_axis(rankdata, 0, X)
    return (R - R.mean(0)) / R.std(0)


def main():
    args = parse_args()
    names = args.collection or ["scie", "emt", "gavish_tnbc"]
    colls = {"scie": SC.SCIE, "emt": SC.EMT, "gavish_tnbc": SC.GAVISH_TNBC}
    C.banner("the null of Route A's own statistic")

    embed = ad.read_h5ad(C.EMBED_H5AD)
    # Live dimensions only: vanished ones would add their near-zero rho to the null.
    live = ~embed.var["vanished"].astype(bool).to_numpy()
    L = ranked_z(np.asarray(embed.X)[:, live])
    n_cells = L.shape[0]
    print(f"{C.RUN_ID}: {n_cells:,} cells x {L.shape[1]} non-vanished dimensions "
          f"= {2 * L.shape[1]} directions")
    rng = np.random.default_rng(args.seed)

    summary = []
    for name in names:
        coll = colls[name]
        sets = C.read_gmt(C.gmt_path(coll))
        cache = CS.tum_dir() / f"random_signature_scores_{coll.name}_{C.RUN_ID}.csv"
        rnd = SP.build_random_scores(coll, sets, embed.obs_names.astype(str),
                                     embed.obs["cohort"], args.n_random, cache,
                                     overwrite=False, seed=args.seed)
        if rnd.empty:
            print(f"[skip] {name}: no matched random sets")
            continue
        base = sorted({c.split("__rnd")[0] for c in rnd.columns})
        draws = sorted({int(c.split("__rnd")[1]) for c in rnd.columns})
        R = pd.DataFrame((L.T @ ranked_z(rnd.to_numpy())) / n_cells,
                         index=embed.var["title"].to_numpy()[live], columns=rnd.columns)

        vals = []
        for _ in range(N_PSEUDO):
            cols = [f"{b}__rnd{rng.choice(draws)}" for b in base]
            sub = R[cols].to_numpy()
            vals.append(sub.max(axis=1))          # the '+' side
            vals.append((-sub).max(axis=1))       # the '-' side
        null = np.concatenate(vals)

        real = pd.read_csv(C.table_path("convergence", coll), comment="#")["A_rho"]
        sizes = [len(v) for v in sets.values()]
        row = {"collection": name, "n_gene_sets": len(base),
               "median_genes": int(np.median(sizes)), "n_random_per_set": len(draws),
               "real_A_rho_p50": float(real.median()), "real_A_rho_max": float(real.max())}
        row.update({f"null_p{p}": float(np.percentile(null, p)) for p in PCTS})
        row["null_max"] = float(null.max())
        for bar in BARS:
            row[f"pass_{bar:g}"] = float((null >= bar).mean())
            row[f"real_{bar:g}"] = int((real >= bar).sum())
        summary.append(row)

        print(f"\n{name}: {len(base)} gene sets (median {int(np.median(sizes))} mapped genes) "
              f"x {len(draws)} matched random each")
        print("  null    " + "  ".join(f"p{p:g} {np.percentile(null, p):.3f}" for p in PCTS)
              + f"  max {null.max():.3f}")
        print(f"  real    p50 {real.median():.3f}  max {real.max():.3f}")
        print("  bar      real directions      by chance")
        for bar in BARS:
            fp = float((null >= bar).mean())
            print(f"  {bar:.3f}    {int((real >= bar).sum()):3d} / {len(real)}"
                  f"            {fp * 100:5.2f}%  (~{fp * len(real):.1f} of {len(real)})")

    if summary:
        df = pd.DataFrame(summary).set_index("collection")
        for name in df.index:
            C.write_table(df.loc[[name]], "route_a_null", colls[name])
        C.banner("side by side")
        cols = ["n_gene_sets", "median_genes", "real_A_rho_p50", "real_A_rho_max",
                "null_p50", "null_p95", "null_p97.5", "null_p99", "null_max"]
        print(df[cols].to_string(float_format="%.3f"))


if __name__ == "__main__":
    main()
