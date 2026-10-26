"""03_2: the cells of this phase and the names of its files.

Malignant cells, 2,000 HVGs without MT- genes, 64 latent dimensions: run id drvi_tum_64_nomt.
Every path and threshold of the phase is defined here, under $DATA_DIR/03_tum/.
"""

from __future__ import annotations

import os
from pathlib import Path

MALIGNANT_LABEL = "malignant"

# Immune labels of the CellTypist model, as in 01_4/fraction_reassignment.py.
IMMUNE_CELL_TYPES = {
    "CD4-Tem", "CD4-Th", "CD4-Th-like", "CD4-Treg", "CD4-activated", "CD4-naive",
    "CD8-Tem", "CD8-Trm", "CD8-activated", "GD", "NK", "NK-ILCs", "NKT", "T_prol",
    "Macro-IFN", "Macro-lipo", "Macro-m1", "Macro-m1-CCL", "Macro-m2", "Macro-m2-CXCL",
    "Mast", "Mono-classical", "Mono-non-classical", "Neutrophil",
    "cDC1", "cDC2", "mDC", "pDC", "mye-prol",
    "b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG",
}


def reassign_fraction(labels):
    """`fraction` recomputed from a post-CNV `cell_type`, with 01_4's rule."""
    import pandas as pd
    return pd.Categorical(
        ["imm" if v in IMMUNE_CELL_TYPES else "non_imm" for v in labels],
        categories=["imm", "non_imm"],
    )

# `cell_type` is constant ('malignant'), so the pre-CNV label is kept as the only real grouping.
LABEL_KEY = "cell_type"
PRIOR_LABEL_KEY = "cell_type_01_4" # the 01_4 CellTypist label
STATUS_KEY = "cnv_status"
COMPARTMENT_KEY = "compartment"
BATCH_KEY = "cohort"
TREATMENT_KEY = "treatment"

# DRVI needs enough cells per cohort (its batch key) to estimate the batch parameters.
MIN_CELLS_PER_COHORT = 100

# Cohorts missing a timepoint are kept, not dropped for one empty biopsy.
DROP_INCOMPLETE_COHORTS = False
REQUIRED_TREATMENTS = ("BASE", "PD1", "RTPD1")

MIN_GENES = 100 # cell filter
MAX_PCT_MT = 10 # cell filter
MIN_CELLS = 3 # gene filter
N_HVGS = 2000
SEED = 0

PREFIX = "shiao_tum" # file prefix of every object of the phase
COMPARTMENT = "tum" # obs['compartment'] and the run id segment
HVG_TAG = "_nomt" # the HVG panel without the MT- genes (reduce_data_tum.py)
N_LATENT = 64 # DRVI latent dimensions

# Leiden on the DRVI latent space, computed in 03_3; `drvi_leiden` is the resolution 03_4 - 03_9 read.
LEIDEN_KEY = "drvi_leiden"
LEIDEN_RESOLUTIONS = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
LEIDEN_RESOLUTION = 1.5 # the resolution 03_4 - 03_9 read, as drvi_leiden
LEIDEN_K = 15 # neighbours, as the UMAP graph of 03_3
LEIDEN_MIN_CELLS = 100 # smaller clusters are merged into their neighbours (isolated ones are kept)


def status_csv() -> Path:
    """The 03_1 malignant call."""
    return tum_dir() / "cnv_status.csv"


def annot_csv() -> Path:
    """The 03_1 CellTypist re-run."""
    return tum_dir() / "cell_annotation.csv"


def prefix() -> str:
    return PREFIX


def compartment() -> str:
    return COMPARTMENT


def data_dir() -> Path:
    """$DATA_DIR, never defaulted: the scripts abort if it is unset."""
    return Path(os.environ["DATA_DIR"]).expanduser().resolve()


def tum_dir() -> Path:
    """$DATA_DIR/03_tum, shared with 03_1."""
    d = data_dir() / "03_tum"
    d.mkdir(parents=True, exist_ok=True)
    return d


def path(suffix: str) -> Path:
    """`path('_norm.h5ad')` -> $DATA_DIR/03_tum/shiao_tum_norm.h5ad."""
    return tum_dir() / f"{PREFIX}{suffix}"


def hvg_tag() -> str:
    return HVG_TAG


def hvg_path(suffix: str) -> Path:
    """`hvg_path('.h5ad')` -> $DATA_DIR/03_tum/shiao_tum_hvg_2k_nomt.h5ad."""
    return path(f"_hvg_2k{HVG_TAG}{suffix}")


def n_latent() -> int:
    return N_LATENT


def leiden_key(res: float) -> str:
    """`leiden_key(0.5)` -> 'drvi_leiden_0.5'."""
    return f"{LEIDEN_KEY}_{res}"


def run_id() -> str:
    """drvi_tum_64_nomt: the id every 03_3 - 03_9 output is named after."""
    return f"drvi_{COMPARTMENT}_{N_LATENT}{HVG_TAG}"


def banner(step: str) -> None:
    """The same header in every log."""
    print(f"{step}  |  run {run_id()}", flush=True)
    print(f"DATA_DIR : {data_dir()}", flush=True)
    print(flush=True)
