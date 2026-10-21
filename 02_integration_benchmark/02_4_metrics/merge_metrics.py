#!/usr/bin/env python
"""
02_4 metrics - merge the per-run CSVs into the single table the plotting reads

One row per <root>/<task>/metrics/<scaling>/<hvg>/<method>_<type>.csv (path
relative to --root, as plotSingleTaskRNA.R expects), one column per metric.
Only SCORED_METRICS are kept: scib also returns rows for metrics not computed here.
PCR batch is blanked on scaled runs, also in CSVs scored before metrics.py skipped it.

Usage:
  python3 merge_metrics.py -o <merged.csv> -r <root/> -i <csv> [<csv> ...]
  python3 merge_metrics.py -o <merged.csv> -r <root/> --glob '<root>/**/hvg/*.csv'   # or let it glob the tree itself
"""

import argparse
import glob as globmod
from functools import reduce

import pandas as pd

# Same as metrics_shared.ALL_METRICS, copied to avoid importing scib and R here.
SCORED_METRICS = (
    "NMI_cluster/label", "ARI_cluster/label", "ASW_label", "ASW_label/batch",
    "PCR_batch", "cell_cycle_conservation", "isolated_label_F1",
    "isolated_label_silhouette", "graph_conn", "iLISI", "cLISI",
    "hvg_overlap",
)


def parse_args():
    p = argparse.ArgumentParser(description="Merge per-run metric CSVs into one table")
    p.add_argument("-o", "--output", required=True, help="merged CSV")
    p.add_argument("-r", "--root", required=True,
                   help="root the column names are made relative to (e.g. .../02_metrics/)")
    p.add_argument("-i", "--input", nargs="*", default=[], help="the per-run CSVs")
    p.add_argument("--glob", default=None, help="glob pattern to collect the CSVs instead of listing them")
    return p.parse_args()


def main():
    args = parse_args()

    files = list(args.input)
    if args.glob:
        files += sorted(globmod.glob(args.glob, recursive=True))
    files = sorted(set(files))
    if not files:
        raise SystemExit("no input CSVs given (use -i or --glob)")

    root = args.root if args.root.endswith("/") else args.root + "/"

    res_list = []
    for file in files:
        clean = file.replace(root, "").replace(".csv", "").lstrip("/")
        res = pd.read_csv(file, index_col=0)
        # Warn rather than silently lose a metric.
        absent = [m for m in SCORED_METRICS if m not in res.index]
        if absent:
            print(f"[warn] {file}: no row for {', '.join(absent)}", flush=True)
        res = res.loc[[m for m in SCORED_METRICS if m in res.index]]
        # 0 by construction on a scaled reference (see metrics_shared.metric_flags).
        if "/scaled/" in f"/{clean}/" and "PCR_batch" in res.index:
            res.loc["PCR_batch"] = float("nan")
        res.rename(columns={res.columns[0]: clean}, inplace=True)
        res_list.append(res)

    results = reduce(
        lambda left, right: pd.merge(left, right, left_index=True, right_index=True),
        res_list,
    )
    results = results.T
    results.to_csv(args.output)
    print(f"merged {len(files)} CSV(s) -> {args.output}  ({results.shape[0]} runs x {results.shape[1]} metrics)",
          flush=True)


if __name__ == "__main__":
    main()
