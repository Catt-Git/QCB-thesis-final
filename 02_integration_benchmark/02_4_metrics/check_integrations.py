#!/usr/bin/env python
"""
02_4 metrics - pre-flight check on an integrated object

Checks in seconds what would otherwise fail hours into metrics.py (no scib):
  * same cells as the reference, in the same order;
  * batch_key and label_key on the reference;
  * the representation each type needs, finite:
      full  -> .X
      embed -> obsm['X_emb']
      knn   -> obsp connectivities + distances

Usage:
  python3 check_integrations.py -i <integrated.h5ad> -u <reference.h5ad> \
      --types full,embed [-b cohort -l cell_type]

Exit code 0 if every check passes, 1 otherwise.
"""

import argparse
import sys

import numpy as np
import anndata as ad
from scipy import sparse


def _finite(matrix) -> bool:
    """True if every entry of a dense or sparse matrix is finite."""
    data = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    return bool(np.isfinite(data).all())


def parse_args():
    p = argparse.ArgumentParser(description="Pre-flight check on an integrated object")
    p.add_argument("-i", "--integrated", required=True)
    p.add_argument("-u", "--uncorrected", required=True, help="reference .h5ad (matching scaling)")
    p.add_argument("--types", required=True,
                   help="comma-separated scib output type(s): knn, embed, full, full,embed")
    p.add_argument("-b", "--batch_key", default="cohort")
    p.add_argument("-l", "--label_key", default="cell_type")
    return p.parse_args()


def main():
    args = parse_args()
    types = [t.strip() for t in args.types.split(",") if t.strip()]
    problems = []
    notes = []

    print(f"[check] reference : {args.uncorrected}", flush=True)
    adata = ad.read_h5ad(args.uncorrected)
    print(f"[check] integrated: {args.integrated}  (types: {', '.join(types)})", flush=True)
    adata_int = ad.read_h5ad(args.integrated)

    # cell count and order 
    if adata.n_obs != adata_int.n_obs:
        problems.append(f"cell count: reference {adata.n_obs} vs integrated {adata_int.n_obs}")
    else:
        ref_names = adata.obs_names.to_numpy()
        int_names = adata_int.obs_names.to_numpy()
        if np.array_equal(ref_names, int_names):
            print(f"  ok   cell order preserved ({adata.n_obs} cells)", flush=True)
        elif set(ref_names) == set(int_names):
            first = int(np.argmax(ref_names != int_names))
            problems.append(
                f"cell order: same cells but permuted; first mismatch at position "
                f"{first} ({ref_names[first]!r} vs {int_names[first]!r})"
            )
        else:
            missing = len(set(ref_names) - set(int_names))
            extra = len(set(int_names) - set(ref_names))
            problems.append(f"cell identity: {missing} in reference only, {extra} in integrated only")

    # batch / label keys 
    for key in (args.batch_key, args.label_key):
        if key not in adata.obs:
            problems.append(f"reference is missing obs[{key!r}]")
        else:
            print(f"  ok   reference has obs[{key!r}] ({adata.obs[key].nunique()} levels)", flush=True)
        if key not in adata_int.obs:
            notes.append(f"integrated object has no obs[{key!r}] (metrics.py copies it from the reference)")

    # representation per type
    for type_ in types:
        if type_ == "full":
            if adata_int.X is None:
                problems.append("full: .X is missing")
            elif not _finite(adata_int.X):
                problems.append("full: .X contains non-finite values")
            else:
                print("  ok   full: .X present and finite", flush=True)
        elif type_ == "embed":
            if "X_emb" not in adata_int.obsm:
                problems.append("embed: obsm['X_emb'] is missing")
            elif not _finite(adata_int.obsm["X_emb"]):
                problems.append("embed: obsm['X_emb'] contains non-finite values")
            else:
                print(f"  ok   embed: obsm['X_emb'] present and finite "
                      f"({adata_int.obsm['X_emb'].shape[1]} dims)", flush=True)
        elif type_ == "knn":
            missing = [k for k in ("connectivities", "distances") if k not in adata_int.obsp]
            if missing:
                problems.append(f"knn: obsp missing {missing}")
            else:
                print("  ok   knn: obsp connectivities + distances present", flush=True)
        else:
            problems.append(f"unknown output type {type_!r}")

    # verdict 
    for note in notes:
        print(f"  note {note}", flush=True)

    if problems:
        print("\n[check] FAILED:", flush=True)
        for prob in problems:
            print(f"  - {prob}", flush=True)
        return 1

    print("\n[check] PASSED: object is a valid metrics input.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
