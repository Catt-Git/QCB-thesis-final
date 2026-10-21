#!/usr/bin/env python
"""
02_4 metrics - the 12 scib metrics, one job per (run, output type)

Scores one integrated object against its reference and writes a one-column CSV
(<method>_<type>). The preparation before scib lives in ../utils/metrics_shared.py.
NMI/ARI use igraph's leiden (leidenalg takes ~20 h here); --hvgs 0 skips HVG
re-selection, since the input is already HVG-only.

Usage:
  python3 metrics.py -u <reference.h5ad> -i <integrated.h5ad> -o <out.csv> \
      -m harmony --type embed [--scaling unscaled -b cohort -l cell_type --organism human --hvgs 0 -v]

The -u reference must match the integrated object's scaling variant; --scaling scaled
skips PCR batch (0 by construction on a per-batch scaled reference).
"""

import argparse
import os
import sys

# metrics_shared imports scib_compat before scib.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import metrics_shared as shared  # noqa: E402
import scib  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Compute the 12 scib metrics")
    p.add_argument("-u", "--uncorrected", required=True, help="reference .h5ad (matching scaling)")
    p.add_argument("-i", "--integrated", required=True, help="integrated .h5ad")
    p.add_argument("-o", "--output", required=True, help="output CSV")
    p.add_argument("-m", "--method", required=True, help="method name, for the column <method>_<type>")
    p.add_argument("--type", required=True, choices=shared.RESULT_TYPES, help="scib output type")
    p.add_argument("--scaling", choices=("unscaled", "scaled"), default="unscaled",
                   help="scaling variant of the run; scaled skips PCR batch")
    p.add_argument("-b", "--batch_key", default="cohort")
    p.add_argument("-l", "--label_key", default="cell_type")
    p.add_argument("--organism", default="human")
    p.add_argument("--hvgs", type=int, default=0,
                   help="HVGs for scib; 0 means no re-selection (the reference is already HVG-restricted)")
    p.add_argument("--cluster-flavor", choices=("igraph", "leidenalg"), default="igraph",
                   help="leiden implementation for the NMI/ARI clustering. igraph is orders of "
                        "magnitude faster and is scanpy's recommended backend; leidenalg is scib's "
                        "own default. Must be the same for every run in the benchmark.")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    type_ = args.type
    setup = f"{args.method}_{type_}"
    n_hvgs = args.hvgs if args.hvgs > 0 else None

    if args.verbose:
        print("=" * 70, flush=True)
        print(f"METRICS  {setup}", flush=True)
        print(f"  reference : {args.uncorrected}", flush=True)
        print(f"  integrated: {args.integrated}", flush=True)
        print(f"  batch/label: {args.batch_key} / {args.label_key}", flush=True)
        print(f"  organism  : {args.organism}   hvgs: {n_hvgs}   scaling: {args.scaling}", flush=True)
        print(f"  clustering: leiden ({args.cluster_flavor})", flush=True)
        print("=" * 70, flush=True)

    # Bind the leiden flavor before scib looks up sc.tl.leiden.
    if args.cluster_flavor == "igraph":
        shared.use_igraph_leiden(verbose=args.verbose)

    adata, adata_int = shared.load_and_align(
        args.uncorrected, args.integrated, args.batch_key, args.label_key, verbose=args.verbose
    )
    shared.ensure_reference_reduced(adata, verbose=args.verbose)
    shared.reduce_integrated(adata_int, type_, n_hvgs=n_hvgs, verbose=args.verbose)

    flags = shared.metric_flags(type_, args.scaling)
    embed = shared.embed_key(type_)

    # The optimal-resolution clustering NMI/ARI need is written next to the CSV.
    file_stump = os.path.splitext(args.output)[0]
    cluster_nmi = f"{file_stump}_nmi.txt"

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    print(f"computing metrics for {setup}", flush=True)
    results = scib.me.metrics(
        adata,
        adata_int,
        batch_key=args.batch_key,
        label_key=args.label_key,
        embed=embed,
        type_=type_,
        organism=args.organism,
        cluster_nmi=cluster_nmi,
        nmi_method="arithmetic",
        nmi_dir=None,
        n_isolated=None,
        verbose=args.verbose,
        **flags,
    )
    results.rename(columns={results.columns[0]: setup}, inplace=True)
    results.to_csv(args.output)

    if args.verbose:
        print(results, flush=True)
    print(f"done -> {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
