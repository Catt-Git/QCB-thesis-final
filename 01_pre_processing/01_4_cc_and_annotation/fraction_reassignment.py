"""
01_4 fraction_reassignment: recode `fraction` to a biological immune label from CellTypist

fraction stored the CD45 sort, which duplicates dataset_origin. It is recoded to 'imm' /
'non_imm' from the CellTypist cell_type, because many cells disagree with their sort
(~10.7% of CD45+ get a non-immune label, ~4.4% of CD45- an immune one, see 01_6).
The sort stays in dataset_origin.

The CD45- fraction is missing for about half of the patients (18/34, see 01_1), while the
processed data uploaded with the paper has non-immune cells for every patient. The authors
probably did a similar reassignment, but this is not stated in the paper :(

The lineage sets below match the Cells_Adult_Breast.pkl model (58 labels).
'Lymph-*' are lymphatic endothelial cells, not lymphocytes, so they are non-immune.

Input : $DATA_DIR/all_samples_combined_scrublet_norm_cc_annotated.h5ad
        Output of celltypist_annotation.py. Must contain .obs['cell_type'] and .obs['fraction'].
Output: same file (in-place). Only .obs['fraction'] changes (CD45+/CD45- -> imm/non_imm);
        everything else, including .obs['dataset_origin'], is untouched.

Usage:
  export DATA_DIR=<data root>
  python3 fraction_reassignment.py
"""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm_cc_annotated.h5ad")
OUT_PATH = IN_PATH  # in-place

# CellTypist labels (Cells_Adult_Breast.pkl, 58 total) grouped by lineage.
IMMUNE_CELL_TYPES = {
    # T / NK
    "CD4-Tem", "CD4-Th", "CD4-Th-like", "CD4-Treg", "CD4-activated", "CD4-naive",
    "CD8-Tem", "CD8-Trm", "CD8-activated", "GD", "NK", "NK-ILCs", "NKT", "T_prol",
    # Myeloid
    "Macro-IFN", "Macro-lipo", "Macro-m1", "Macro-m1-CCL", "Macro-m2", "Macro-m2-CXCL",
    "Mast", "Mono-classical", "Mono-non-classical", "Neutrophil",
    "cDC1", "cDC2", "mDC", "pDC", "mye-prol",
    # B / Plasma
    "b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG",
}
NON_IMMUNE_CELL_TYPES = {
    # Epithelial (luminal hormone-responsive / luminal secretory / basal)
    "LummHR-SCGB", "LummHR-active", "LummHR-major",
    "Lumsec-HLA", "Lumsec-KIT", "Lumsec-basal", "Lumsec-lac", "Lumsec-major",
    "Lumsec-myo", "Lumsec-prol", "basal",
    # Stromal / vascular
    "Fibro-SFRP4", "Fibro-major", "Fibro-matrix", "Fibro-prematrix",
    "Lymph-immune", "Lymph-major", "Lymph-valve1", "Lymph-valve2",
    "Vas-arterial", "Vas-capillary", "Vas-venous", "pericytes", "vsmc",
}

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

assert "cell_type" in adata.obs, "Expected CellTypist labels in .obs['cell_type']"
assert "fraction" in adata.obs, "Expected .obs['fraction']"

# Every cell_type must be in one of the two sets, otherwise unknown labels end up in 'non_imm'.
observed = set(adata.obs["cell_type"].astype(str).unique())
unclassified = observed - IMMUNE_CELL_TYPES - NON_IMMUNE_CELL_TYPES
assert not unclassified, (
    f"{len(unclassified)} cell_type value(s) not in either immune/non-immune list; "
    f"update the lineage sets before recoding: {sorted(unclassified)}"
)

is_immune_celltype = adata.obs["cell_type"].isin(IMMUNE_CELL_TYPES)

old_fraction = adata.obs["fraction"].copy()
adata.obs["fraction"] = pd.Categorical(
    np.where(is_immune_celltype, "imm", "non_imm"),
    categories=["imm", "non_imm"],
)

print("\nNew `fraction` counts:", flush=True)
print(adata.obs["fraction"].value_counts(), flush=True)
print("\nOld CD45 sort (fraction) vs new CellTypist-based fraction:", flush=True)
print(pd.crosstab(old_fraction, adata.obs["fraction"]), flush=True)
print("\ndataset_origin (unchanged) vs new fraction:", flush=True)
print(pd.crosstab(adata.obs["dataset_origin"], adata.obs["fraction"]), flush=True)

# Post-conditions: fraction fully recoded, technical sort preserved.
assert set(adata.obs["fraction"].unique()) <= {"imm", "non_imm"}
assert adata.obs["fraction"].notna().all(), "recoded fraction has NaNs"
assert "dataset_origin" in adata.obs, "dataset_origin must be preserved"

# Write to a temp file, then replace: an interrupted write does not corrupt the input file.
print("\nSaving (in-place, atomic)...", flush=True)
tmp_path = OUT_PATH + ".tmp"
adata.write_h5ad(tmp_path, compression="gzip")
os.replace(tmp_path, OUT_PATH)
print(f"Wrote {OUT_PATH}", flush=True)
