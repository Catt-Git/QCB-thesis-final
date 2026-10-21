"""
02 utils - write .h5ad files that older anndata versions can still read

Inputs are written with anndata 0.13 / pandas 3, but scgen-py is pinned to
anndata 0.10.8, which cannot read these encodings:

    nullable-string-array   obs/_index, var/_index, var/gene_ids
    nullable-boolean        obs/predicted_doublet
    null                    uns/log1p/base   (written by scanpy)

They are downcast before writing (lossless, still read by anndata 0.13).

Usage:
  from h5ad_compat import write_h5ad_compat
  write_h5ad_compat(adata, path, compression="gzip")   # instead of adata.write_h5ad()

The written file is re-checked at HDF5 level to catch encodings not handled yet.
"""

from __future__ import annotations

import h5py
import numpy as np
import pandas as pd

# Encodings anndata 0.10.x cannot read (nullable-integer not seen, but guarded).
LEGACY_INCOMPATIBLE_ENCODINGS = frozenset(
    {"nullable-string-array", "nullable-boolean", "nullable-integer", "null"}
)


def _sanitize_frame(df: pd.DataFrame, name: str, changed: list[str]) -> None:
    """Downcast index and extension-dtype columns of one of .obs / .var."""
    if df.index.dtype != object:
        df.index = df.index.astype(object)
        changed.append(f"{name}.index -> object")

    for column in df.columns:
        dtype = df[column].dtype

        # Categoricals are written as `categorical`, which every version reads.
        if isinstance(dtype, pd.CategoricalDtype):
            continue

        if isinstance(dtype, pd.StringDtype):
            df[column] = df[column].astype(object)
            changed.append(f"{name}[{column!r}] string -> object")

        elif pd.api.types.is_extension_array_dtype(dtype):
            # Nullable dtypes: downcast only if there are no missing values.
            n_missing = int(df[column].isna().sum())
            if n_missing:
                raise ValueError(
                    f"{name}[{column!r}] is {dtype} with {n_missing} missing "
                    "values; downcasting to a numpy dtype would have to invent a "
                    "value for them. Decide explicitly what those cells mean "
                    "before writing."
                )
            df[column] = df[column].to_numpy()
            changed.append(f"{name}[{column!r}] {dtype} -> numpy")


def _drop_none(mapping, path: str, changed: list[str]) -> None:
    """Remove None values from .uns, recursively (e.g. scanpy's uns['log1p']['base'])."""
    for key in list(mapping):
        value = mapping[key]
        if isinstance(value, dict):
            _drop_none(value, f"{path}[{key!r}]", changed)
        elif value is None:
            del mapping[key]
            changed.append(f"{path}[{key!r}] None -> removed")


def sanitize_for_legacy_readers(adata) -> list[str]:
    """Rewrite the dtypes that only anndata >=0.11 can read. Mutates in place.

    :param adata: the object about to be written
    :return: human-readable list of the changes made, possibly empty
    """
    changed: list[str] = []
    _sanitize_frame(adata.obs, "obs", changed)
    _sanitize_frame(adata.var, "var", changed)
    _drop_none(adata.uns, "uns", changed)
    return changed


def find_incompatible_encodings(path) -> list[tuple[str, str]]:
    """List every dataset in a written .h5ad an old anndata could not read.

    :param path: path to an .h5ad file
    :return: list of (dataset path, encoding) pairs, empty if the file is clean
    """
    hits: list[tuple[str, str]] = []

    def visit(name, obj):
        encoding = obj.attrs.get("encoding-type")
        if encoding in LEGACY_INCOMPATIBLE_ENCODINGS:
            hits.append((name, encoding))

    with h5py.File(path, "r") as handle:
        handle.visititems(visit)
        # visititems does not visit the root, which carries the AnnData encoding.
        root = handle.attrs.get("encoding-type")
        if root in LEGACY_INCOMPATIBLE_ENCODINGS:
            hits.append(("/", root))

    return hits


def assert_legacy_readable(path) -> None:
    """Fail if the written file still contains an encoding old readers reject."""
    hits = find_incompatible_encodings(path)
    if hits:
        listing = "\n".join(f"    {encoding:<22} {name}" for name, encoding in hits)
        raise AssertionError(
            f"{path} contains encodings that anndata 0.10 (the scgen-py "
            f"environment) cannot read:\n{listing}\n"
            "Extend sanitize_for_legacy_readers() in utils/h5ad_compat.py to "
            "cover them."
        )


def write_h5ad_compat(adata, path, verbose: bool = True, **kwargs) -> list[str]:
    """Sanitise, write, and verify what landed on disk.

    Drop-in replacement for ``adata.write_h5ad(path, **kwargs)``.

    :param adata: object to write; modified in place by the sanitiser
    :param path: destination .h5ad
    :param verbose: print the changes made
    :param kwargs: forwarded to ``write_h5ad`` (e.g. compression)
    :return: the list of changes made
    """
    changed = sanitize_for_legacy_readers(adata)
    if verbose:
        if changed:
            print(f"[h5ad_compat] {len(changed)} dtype adjustments before writing:", flush=True)
            for entry in changed:
                print(f"    {entry}", flush=True)
        else:
            print("[h5ad_compat] nothing to adjust", flush=True)

    adata.write_h5ad(path, **kwargs)
    assert_legacy_readable(path)

    if verbose:
        print(f"[h5ad_compat] wrote {path}, readable by anndata >=0.7", flush=True)
    return changed
