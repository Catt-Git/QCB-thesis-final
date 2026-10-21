#!/usr/bin/env python
"""
02_2 integration: scGen, in its own environment.

Separate script because scgen-py (scvi-tools 0.19, torch 2.0) cannot coexist
with the main stack and has no scib: SCGEN is called directly with scib's call
sequence (setup_anndata -> SCGEN -> train -> batch_removal). Output type: full.

Added on top of scGen: a model checkpoint in --model-dir, and a chunked decoder
pass in batch_removal (619k cells do not fit on a 4 GB GPU).

Usage:
  conda activate scgen-py
  export DATA_DIR=<data root>
  python3 run_scgen.py \
      -i $DATA_DIR/shiao_hvg_2k.h5ad \
      -o $DATA_DIR/02_integration/scgen_unscaled.h5ad
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import torch
import anndata as ad
from scgen import SCGEN

from model_paths import (default_model_dir, has_saved_model, model_prefix,
                         saved_model)


def parse_args():
    p = argparse.ArgumentParser(description="Run scGen integration")
    p.add_argument("-i", "--input", required=True, help="prepared .h5ad input")
    p.add_argument("-o", "--output", required=True, help="integrated .h5ad output")
    p.add_argument("-b", "--batch-key", default="cohort")
    p.add_argument("-l", "--label-key", default="cell_type")
    p.add_argument("--epochs", type=int, default=100,
                   help="training epochs [default %(default)s, scib's value]")
    p.add_argument("--model-dir", default=None,
                   help="directory the trained model is saved to and reloaded "
                        "from, as <run_id>_model.pt "
                        "[default: <DATA_DIR>/02_scgen]")
    p.add_argument("--retrain", action="store_true",
                   help="train even if --model-dir holds a saved model")
    p.add_argument("--decode-chunk", type=int, default=16384,
                   help="cells per decoder forward pass in batch_removal "
                        "[default %(default)s]")
    return p.parse_args()


def _chunk_generative(model, chunk_size):
    """Make the decoder pass inside batch_removal run in chunks.

    batch_removal() decodes all cells at once (OOM on a 4 GB GPU). The decoder is
    row-wise in eval mode, so chunking gives the same result with less memory.
    """
    module = model.module
    original = module.generative          # bound, keeps scvi's arg-transfer decorator

    def generative(z, *args, **kwargs):
        if z.dim() != 2 or z.shape[0] <= chunk_size:
            return original(z, *args, **kwargs)
        print(f"[decode] {z.shape[0]:,} cells in chunks of {chunk_size:,}", flush=True)
        out = []
        for start in range(0, z.shape[0], chunk_size):
            px = original(z[start:start + chunk_size], *args, **kwargs)["px"]
            out.append(px.detach().to("cpu"))
        return {"px": torch.cat(out)}

    module.generative = generative


def _fit_or_load(adata, args):
    """Return a trained SCGEN, reusing a saved checkpoint when there is one."""
    SCGEN.setup_anndata(adata, batch_key=args.batch_key, labels_key=args.label_key)
    prefix = model_prefix(args.output)
    checkpoint = saved_model(args.model_dir, prefix)

    if not args.retrain and has_saved_model(args.model_dir, prefix):
        print(f"[load] trained model from {checkpoint} (--retrain to ignore)", flush=True)
        return SCGEN.load(args.model_dir, adata=adata, prefix=prefix)

    # Same parametrisation as scib.integration.scgen (from the scGen tutorial).
    print(f"[run] scGen (batch_key={args.batch_key}, labels_key={args.label_key}, "
          f"epochs={args.epochs})", flush=True)
    model = SCGEN(adata)
    model.train(
        max_epochs=args.epochs,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25,
    )
    print(f"[save] trained model -> {checkpoint}", flush=True)
    model.save(args.model_dir, prefix=prefix, overwrite=True)
    return model


def main():
    args = parse_args()
    if args.model_dir is None:
        args.model_dir = default_model_dir(args.output, "scgen")

    print(f"[read] {args.input}", flush=True)
    adata = ad.read_h5ad(args.input)
    print(f"[read] {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

    for key in (args.batch_key, args.label_key):
        assert key in adata.obs, f"missing obs column {key!r}"

    obs_names_before = adata.obs_names.to_numpy().copy()

    model = _fit_or_load(adata, args)

    _chunk_generative(model, args.decode_chunk)
    print("[run] batch_removal", flush=True)
    corrected = model.batch_removal()

    # batch_removal returns a new object; make sure it matches cell for cell.
    assert corrected.n_obs == adata.n_obs, "scGen changed the cell count"
    assert np.array_equal(corrected.obs_names.to_numpy(), obs_names_before), (
        "scGen changed the cell order; the metrics would be computed on "
        "mismatched cells"
    )
    assert np.isfinite(corrected.X).all(), "scGen produced non-finite values"

    # Same output shape as run_integration.py: no raw counts.
    if "counts" in corrected.layers:
        del corrected.layers["counts"]
    # Drop scGen's own latent ('latent' or 'corrected_latent') so it is not taken for X_emb.
    for latent_key in ("latent", "corrected_latent"):
        corrected.obsm.pop(latent_key, None)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    print(f"[write] {args.output}", flush=True)
    corrected.write_h5ad(args.output, compression="gzip")

    size_gb = os.path.getsize(args.output) / 1024 ** 3
    print(f"[done] {corrected.n_obs:,} x {corrected.n_vars:,}, {size_gb:.2f} GB", flush=True)


if __name__ == "__main__":
    main()
