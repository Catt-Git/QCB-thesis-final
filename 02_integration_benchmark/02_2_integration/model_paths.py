"""02_2 integration: checkpoint paths for the trained methods (scVI, scANVI, scGen).

Training is the long part, so the model is saved and a crash afterwards does not
redo it. Each method gets its own flat directory, separate from 02_integration/:

    $DATA_DIR/02_scvi/scvi_unscaled_model.pt
    $DATA_DIR/02_scanvi/scanvi_unscaled_scvi_model.pt
    $DATA_DIR/02_scanvi/scanvi_unscaled_scanvi_model.pt
    $DATA_DIR/02_scgen/scgen_unscaled_model.pt
    $DATA_DIR/02_scgen/scgen_scaled_model.pt

Callers pass --model-dir; the run id becomes the file prefix. This module is the
fallback for running by hand. Delete a checkpoint to retrain.
"""

from __future__ import annotations

import os


def run_id_of(output_path):
    """Grid run id of an integration output: <run_id>.h5ad -> <run_id>."""
    return os.path.splitext(os.path.basename(os.path.abspath(output_path)))[0]


def default_model_dir(output_path, method):
    """<DATA_DIR>/02_<method>, with DATA_DIR the output's grandparent."""
    data_dir = os.path.dirname(os.path.dirname(os.path.abspath(output_path)))
    return os.path.join(data_dir, f"02_{method}")


def model_prefix(output_path):
    """File-name prefix for this run's checkpoint(s) inside the model directory."""
    return f"{run_id_of(output_path)}_"


def saved_model(model_dir, prefix):
    """Path scvi-tools' save(prefix=...) writes, or None if there is no directory."""
    if model_dir is None:
        return None
    return os.path.join(model_dir, f"{prefix}model.pt")


def has_saved_model(model_dir, prefix):
    """True if that checkpoint is already on disk."""
    path = saved_model(model_dir, prefix)
    return path is not None and os.path.isfile(path)
