#!/usr/bin/env python
"""
02_2 integration: dispatcher for the Python methods.

Runs one method (bbknn, harmony, scvi, scanvi, scanorama; bodies in
integration_methods.py) on one prepared input. Checks the cell order is unchanged,
drops layers['counts'], and writes the .h5ad plus, for embed outputs, --emb-out .npy.
scVI/scANVI checkpoint to --model-dir (see model_paths.py).

Usage:
  export DATA_DIR=<data root>
  python3 run_integration.py -m harmony \
      -i $DATA_DIR/shiao_hvg_2k.h5ad \
      -o $DATA_DIR/02_integration/harmony_unscaled.h5ad \
      --emb-out $DATA_DIR/02_embeddings/harmony_unscaled.npy
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

# scib_compat must be imported before scib (pulled in by integration_methods).
UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import scib_compat  # noqa: F401,E402

import scanpy as sc  # noqa: E402
import integration_methods as im  # noqa: E402
from model_paths import default_model_dir, model_prefix  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Run one Python integration method")
    p.add_argument("-m", "--method", required=True, choices=sorted(im.METHODS))
    p.add_argument("-i", "--input", required=True, help="prepared .h5ad input")
    p.add_argument("-o", "--output", required=True, help="integrated .h5ad output")
    p.add_argument("-b", "--batch-key", default="cohort")
    p.add_argument("-l", "--label-key", default="cell_type")
    p.add_argument("--emb-out", default=None,
                   help="path to write obsm['X_emb'] as .npy (embed outputs)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--model-dir", default=None,
                   help="directory scVI/scANVI save and reload their trained "
                        "model in, as <run_id>_model.pt; ignored by the other "
                        "methods [default: <DATA_DIR>/02_<method>]")
    p.add_argument("--retrain", action="store_true",
                   help="train even if --model-dir holds a saved model")
    return p.parse_args()


def main():
    args = parse_args()
    method = args.method
    if args.model_dir is None:
        args.model_dir = default_model_dir(args.output, method)

    print(f"[read] {args.input}", flush=True)
    adata = sc.read_h5ad(args.input)
    print(f"[read] {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

    if method in ("scvi", "scanvi"):
        assert "counts" in adata.layers, (
            f"{method} reads raw counts from layers['counts'], which is missing. "
            "Its input must be the unscaled object that preserves that layer."
        )

    obs_names_before = adata.obs_names.to_numpy().copy()

    print(f"[run] {method} (batch_key={args.batch_key}, "
          f"label_key={args.label_key}, seed={args.seed})", flush=True)
    fn = im.METHODS[method]
    integrated = fn(adata, args.batch_key, label_key=args.label_key, seed=args.seed,
                    model_dir=args.model_dir, model_prefix=model_prefix(args.output),
                    retrain=args.retrain)

    # Cell order is load-bearing for the metrics; verify it survived.
    assert np.array_equal(integrated.obs_names.to_numpy(), obs_names_before), (
        f"{method} changed the cell order; the metrics would be computed on "
        "mismatched cells"
    )

    # Report what the object now carries, and sanity-check it against the type.
    has_emb = "X_emb" in integrated.obsm
    has_graph = "connectivities" in integrated.obsp
    print(f"[out] obsm['X_emb']={has_emb}  graph={has_graph}  "
          f"obsm={list(integrated.obsm)}", flush=True)
    if method == "bbknn":
        assert has_graph, "bbknn produced no neighbour graph"
    elif method in ("harmony", "scvi", "scanvi", "scanorama"):
        assert has_emb, f"{method} produced no obsm['X_emb']"

    # The metrics never read raw counts; dropping them keeps the transfer light.
    if "counts" in integrated.layers:
        del integrated.layers["counts"]

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    print(f"[write] {args.output}", flush=True)
    integrated.write_h5ad(args.output, compression="gzip")

    if args.emb_out and has_emb:
        os.makedirs(os.path.dirname(os.path.abspath(args.emb_out)), exist_ok=True)
        np.save(args.emb_out, np.asarray(integrated.obsm["X_emb"], dtype=np.float32))
        print(f"[write] {args.emb_out} ({integrated.obsm['X_emb'].shape})", flush=True)

    size_gb = os.path.getsize(args.output) / 1024 ** 3
    print(f"[done] {integrated.n_obs:,} x {integrated.n_vars:,}, {size_gb:.2f} GB", flush=True)


if __name__ == "__main__":
    main()
