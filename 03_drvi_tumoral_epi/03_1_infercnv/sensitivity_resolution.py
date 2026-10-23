#!/usr/bin/env python3
"""
03_1 sensitivity_resolution: per-patient verdict at each leiden resolution of the sweep

Applies the chain of call_malignant.ipynb to each run: null position of each cell
against its cohort's stromal cells, median per subcluster, subcluster CNV-positive if
the median >= MAL_MED, cohort CNV-positive if >= GATE_FRAC of its epithelium is in such
subclusters. Also reports the share of epithelium in subclusters >= MIN_SUB, which is what
picks the operating resolution. Exits 1 if any cohort's verdict changes.

Input : $DATA_DIR/03_tum/sensitivity/resolution/<cohort>_res<r>.csv
        Output of sensitivity_resolution.sh.
Output: ../tables/03_1_infercnv/sensitivity/resolution_sweep.csv

Usage:
  export DATA_DIR=<data root>
  python3 sensitivity_resolution.py
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np
import pandas as pd

MAL_MED = 0.90 # subcluster median nullpos >= this -> CNV-positive
GATE_FRAC = 0.10 # cohort CNV-positive if this fraction of its epithelium is in such subclusters
MIN_SUB = 20 # reported, not applied


def main() -> int:
    data_dir = os.path.abspath(os.path.expanduser(os.environ["DATA_DIR"]))
    in_dir = os.path.join(data_dir, "03_tum", "sensitivity", "resolution")
    here = os.path.dirname(os.path.abspath(__file__))
    phase = os.path.dirname(here)
    tab_dir = os.path.join(phase, "tables", "03_1_infercnv", "sensitivity")
    os.makedirs(tab_dir, exist_ok=True)

    paths = sorted(glob.glob(os.path.join(in_dir, "*_res*.csv")))
    assert paths, f"no sweep output in {in_dir}; run ./sensitivity_resolution.sh first"

    rows = []
    for p in paths:
        m = re.match(r"(\w+)_res([\d.]+)\.csv", os.path.basename(p))
        cohort, res = m.group(1), float(m.group(2))
        d = pd.read_csv(p)
        stromal = d[d["group"] == "stromal"]
        for ax in ("cnv_score", "cnv_corr"):
            d["p_" + ax] = (np.searchsorted(np.sort(stromal[ax].to_numpy()),
                                            d[ax].to_numpy()) / len(stromal))
        d["nullpos"] = np.minimum(d["p_cnv_score"], d["p_cnv_corr"])

        epi = d[(d["group"] == "epi") & d["subcluster"].notna()]
        g = epi.groupby("subcluster")
        sc = pd.DataFrame({"n": g.size(), "median_nullpos": g["nullpos"].median()})
        sc["positive"] = sc["median_nullpos"] >= MAL_MED
        frac = epi["subcluster"].isin(sc.index[sc["positive"]]).mean()
        big = sc.index[sc["n"] >= MIN_SUB]

        held = d[d["group"] == "immune_heldout"]
        rows.append(dict(
            cohort=cohort, resolution=res, n_epi=len(epi),
            n_subclusters=len(sc), n_subclusters_ge_min=int((sc["n"] >= MIN_SUB).sum()),
            median_subcluster_size=int(sc["n"].median()), min_subcluster_size=int(sc["n"].min()),
            pct_cells_in_subclusters_ge_min=round(100 * epi["subcluster"].isin(big).mean(), 1),
            n_subclusters_positive=int(sc["positive"].sum()),
            frac_epi_positive=round(frac, 4),
            verdict="CNV-positive" if frac >= GATE_FRAC else "CNV-negative",
            median_nullpos_heldout=round(held["nullpos"].median(), 4) if len(held) else np.nan,
            median_nullpos_stromal=round(
                d.loc[d["group"] == "stromal", "nullpos"].median(), 4)))

    r = pd.DataFrame(rows).sort_values(["cohort", "resolution"])
    r.to_csv(os.path.join(tab_dir, "resolution_sweep.csv"), index=False)

    pd.set_option("display.width", 250)
    print(r.to_string(index=False))
    print()
    unstable = []
    for cohort, d in r.groupby("cohort"):
        spread = d["frac_epi_positive"].max() - d["frac_epi_positive"].min()
        stable = d["verdict"].nunique() == 1
        print(f"{cohort}: {'stable' if stable else 'UNSTABLE'} "
              f"{sorted(d['verdict'].unique())} | frac_epi_positive "
              f"{d['frac_epi_positive'].min():.3f}-{d['frac_epi_positive'].max():.3f} "
              f"(spread {spread:.3f})")
        if not stable:
            unstable.append(cohort)

    print(f"\nWrote {tab_dir}/resolution_sweep.csv")
    return 1 if unstable else 0


if __name__ == "__main__":
    sys.exit(main())
