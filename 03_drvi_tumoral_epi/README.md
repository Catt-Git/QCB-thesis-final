# 03_drvi_tumoral_epi

Fourth phase of the thesis: the biological interpretation of the tumor. inferCNV separates the
malignant epithelial cells, DRVI is trained on them alone, and its latent dimensions are read
against three collections of gene signatures (SCIE, EMT, the Gavish pan-cancer metaprograms) 
and across the three treatment timepoints.

Input: `shiao.h5ad` (619,693 × 30,869) from phase 01, read-only. Batch key `cohort`.

## Requirements

> [!IMPORTANT]
> This phase uses **three conda environments**: `benchmark-py-r` for almost everything,
> `infercnv-r` for inferCNV (03_1) and `cytotrace2-py` for CytoTRACE2 (03_5). The drivers
> switch environment where needed; see the README of each step.

- `$DATA_DIR/shiao.h5ad` and `$DATA_DIR/Cells_Adult_Breast.pkl` from phase 01,
  `$DATA_DIR/regev_lab_cell_cycle_genes.txt`, and the signature lists in `$DATA_DIR/signatures/`
  (versioned in `datasets/`).
- DRVI training (03_3) runs locally on a GPU or on the cluster on CPU (`submit_drvi_tum.slurm`).
- Every `.sh` driver (`infercnv_all.sh`, `subsetting_all.sh`, `signatures_all.sh`,
  `cell_first_all.sh`, `factor_first_all.sh`, `convergence_all.sh`, `argmax_all.sh`, the two
  03_1 sensitivity scripts) also takes `--slurm`: the whole driver goes up as one job
  (`utils/submit_driver.slurm`, partition `long`) that runs it with the same arguments and
  resumes the same way. `JOB_ENV`, `JOB_CPUS`, `JOB_MEM` and `JOB_PARTITION` override the defaults.

## Quick start

```bash
cd 03_drvi_tumoral_epi
export DATA_DIR=<data root>
conda env create -f ../environments/infercnv-r.yml # once, for 03_1

cd 03_1_infercnv && bash infercnv_all.sh --threads 12
# run by hand: call_malignant.ipynb
python3 recelltypist_nonmalignant.py

cd ../03_2_subsetting
# run by hand: subset_and_qc.ipynb
bash subsetting_all.sh

cd ../03_3_drvi_run && python3 run_drvi_tum.py # or drvi_tum_64.ipynb, or the checkpoints below

cd ../03_4_signatures && bash signatures_all.sh
cd ../03_5_cytotrace2 && python3 cytotrace2_tum.py  # in the cytotrace2-py env
cd ../03_6_cell_first && bash cell_first_all.sh
cd ../03_7_factor_first && bash factor_first_all.sh
cd ../03_8_convergence && bash convergence_all.sh
cd ../03_8b_argmax && bash argmax_all.sh  # optional

cd ../03_9_treatment && python3 treatment_programs_tum.py
```

From 03_4 to 03_8 each step has its own driver, run in this order (03_7 reads 03_6's row order,
03_8 reads both): each covers all three collections and skips what is already on disk;
`--dry-run` previews, `--force` re-runs, `--slurm` submits it as one job. 03_5 is optional and
runs in its own env (cytotrace2-py pins numpy<2); `cell_first_all.sh` also runs 03_4's
`signature_composition_tum.py`, which needs the 03_6 scores. 03_8b is optional too: it reads
03_8 and nothing downstream reads it.

## Execution order

| # | Folder | Script | What it does |
|---|--------|--------|--------------|
| 1 | `03_1_infercnv` | `infercnv_all.sh`, `call_malignant.ipynb`, `recelltypist_nonmalignant.py` | inferCNV per patient, the malignant call, CellTypist on the other cells |
| 2 | `03_2_subsetting` | `subset_and_qc.ipynb`, `subsetting_all.sh` | the malignant object: QC, scran, cell cycle, HVGs without `MT-` genes |
| 3 | `03_3_drvi_run` | `run_drvi_tum.py` or `drvi_tum_64.ipynb` | DRVI, 64 latent dimensions |
| 4 | `03_4_signatures` | `signatures_all.sh` (`build_signatures_tum.py`, `emt_vs_gavish_tum.py`) | the lists of each collection → `.gmt`, coverage, Jaccard |
| 5 | `03_5_cytotrace2` | `cytotrace2_tum.py` | CytoTRACE2 potency per patient, in its own env |
| 6 | `03_6_cell_first` | `cell_first_all.sh` (`cell_first_tum.py`) | Route A: per-cell signature scores against the dimensions |
| 7 | `03_7_factor_first` | `factor_first_all.sh` (`factor_first_tum.py`) | Route B: enrichment of each dimension's top genes |
| 8 | `03_8_convergence` | `convergence_all.sh` (`convergence_tum.py`) | Route C: where A and B agree, per collection |
| 8 | `03_8_convergence` | `route_a_null_tum.py` | the null the ρ bar is calibrated on |
| 8 | `03_8_convergence` | `leiden_smi_tum.py`, `leiden_states_tum.py` | the DRVI-Leiden clusters: which DR direction owns each (SMI), what they contain (DE, DR enrichment) |
| 8 | `03_8_convergence` | `dr_naming_tum.py`, `dr_umap_tum.py` | which dimensions get a name (with the SMI clusters overlaid), and where |
| 8b | `03_8b_argmax` | `argmax_all.sh` (`argmax_states_tum.py`) | optional: the max(DR) groups, each cell to its top direction, with 03_8's DE and DR enrichment |
| 9 | `03_9_treatment` | `treatment_programs_tum.py` (`treatment_common.py`) | the named programs across the timepoints, paired within the patient |

`utils/signature_common.py` (paths, cutoffs, figure and table writers) and
`utils/sig_collections.py` (the three collections) are shared by 03_4 - 03_9.

## The run

One cell set, one gene panel, one latent space, all fixed in `utils/cell_set.py`:

- **cells:** `cnv_status == 'malignant'` (03_1 call), cohorts with ≥ 100 malignant
  cells: **41,943 cells, 20 cohorts**. Cohorts missing a timepoint are kept.
- **genes:** the 2,000 batch-aware HVGs selected on these cells, with the mitochondrial `MT-`
  genes replaced by the next genes of the same ranking (`reduce_data_tum.py`): they track
  dissociation stress, not a cell state.
- **DRVI:** 64 latent dimensions, `batch_key='cohort'`, no `labels_key`. Run id
  **`drvi_tum_64_nomt`**, in every file name and table header.

Two settings are fixed for 03_4 - 03_9 in `utils/signature_common.py`:

- **vanished dimensions are pruned**: 03_6 drops the dimensions DRVI flagged vanished, 03_7 the
  vanished directions (a dimension can vanish on one side only), 03_8 inherits 03_6's list.
  That is 49 of 64 dimensions for Route A, 97 of 128 directions for Route B (`DR 48+` vanished
  on one side only), 98 rows in the convergence tables.
- **Route A's bar is |ρ| ≥ 0.30**, calibrated on the per-direction null of
  `route_a_null_tum.py`: ~p95.5 on `scie` (null p97.5 = 0.37), ~p96.5 on `gavish_tnbc`, above
  p99 on `emt`, i.e. ~4, ~3 and ~0 of 98 directions pass by chance.

## Design choices

- **cell_type is constant.** After the call every cell is `malignant`, so cells are grouped by
  `drvi_leiden` (Leiden on the live DRVI dimensions, resolution 1.5, clusters under 100 cells
  merged, computed in 03_3) first
  and by `cell_type_01_4`, the pre-CNV CellTypist label, second, as a landmark and never as an
  identity. The Leiden clusters are cut from the DRVI space itself, so a dimension that owns a
  cluster is circular by construction.
- **Route A standardises within `cohort` only.** Inside one malignant compartment the CellTypist
  groups are states of the same tumor; standardising within them (`Lumsec-prol` means
  proliferating) would remove the proliferation axis instead of measuring it.
- **EMT is a continuum.** The `emt` target region is epithelial-high × mesenchymal-high. Being
  aneuploid rules out an actual fibroblast in the object, not ambient RNA.

## The three collections

| collection | slug | question | target region |
|---|---|---|---|
| `scie` | `scie` | stemness × immunogenicity: ten lab lists plus CytoTRACE2 | yes |
| `emt` | `emt` | EMT: collaborator list B, epithelial / hybrid / mesenchymal | yes |
| `gavish` | `gavish_tnbc` | what is each dimension? The Gavish et al. 2023 pan-cancer metaprograms | no |

They share no output: the collection is in the name of every table and figure.

**gavish.** The 41 metaprograms (Gavish et al. 2023, Nature 618:598-606) are in
`$DATA_DIR/signatures/GAVISH_metaprograms/`; `MP1` (G2/M) is not in the MSigDB export and is
taken from the authors' `MP_list.RDS` by `utils/gavish_extraction.py`. By default 22 are scored
(`gavish_tnbc`): the states reported in basal-like / TNBC malignant cells plus `MP7` (in-vitro
stress) to separate stress from dissociation; `--all-metaprograms` scores all 41 under slug
`gavish`. `MP38` (glutathione, proximal-tubule genes) and `MP11` (translation initiation, a
library-complexity pattern) were dropped on a matched random-list null. The set is the tuple
`_GAVISH_TNBC_MPS` in `utils/sig_collections.py`.

Route A and Route B land on the same metaprogram on **25 of the 98 dimension-directions** (27
convergent by family), led by `DR 7- = MP5_STRESS` (ρ 0.814, FDR 6e-40),
`DR 2+ = MP1_CELL_CYCLE_G2_M` (ρ 0.813, FDR 5e-50) and `DR 5+ = MP2_CELL_CYCLE_G1_S` (ρ 0.757).

## Results on `drvi_tum_64_nomt`

| step | result |
|---|---|
| 03_1 | 42,291 malignant cells; 20 cohorts with ≥ 100 of them keep 41,943 |
| 03_3 | 15 of 64 dimensions vanished, 49 live (97 live directions); 34 Leiden clusters at 1.5 |
| 03_6 Route A | \|ρ\| ≥ 0.30 on 16 directions on `scie`, 1 on `emt`, 34 on `gavish_tnbc`; the `scie` target (2,704 cells) is mostly cycling |
| 03_7 Route B | directions with a significant set: 53 on `scie`, 3 on `emt` (mesenchymal only), 88 on `gavish_tnbc` |
| 03_8 convergence | `convergent` rows of 98: 13 on `scie`, 0 on `emt`, 27 on `gavish_tnbc`. No `scie` or `emt` target axis passes both criteria: no stem-like immune-evasive axis, no hybrid E/M axis |
| 03_8 naming | 22 programs named, 38 bars on 30 directions of 28 dimensions; `DR 1+` (`MP12_EMT_1`) is named too but its dimension is a depth axis |
| 03_8 clusters | 7 directions own a Leiden cluster (SMI ≥ 0.7, one to one), none of them named |
| 03_9 | 20 patients, 49 biopsies. A priori: G2/M and G1/S cycle down after RT (`RTPD1_vs_PD1`, 8 of 10 patients, FDR 0.029; only G2/M confirmed at gene level); IFN-γ, MHC-II, type-I IFN and senescence move the expected way, not significant. Exploratory: no test at FDR < 0.05 |

The `CLUSTER_READING` of 03_8 and the `ARGMAX_PLAIN_READING` of 03_8b are drafts (2026-10-27)
still to be checked by eye.

## Checkpoints

To read the reported run without retraining, put these files in `$DATA_DIR/03_tum/` and start
from 03_4:

| File | Written by | What it is |
|---|---|---|
| `model_drvi_tum_64_nomt.pt` | 03_3 | the trained DRVI model |
| `embed_drvi_tum_64_nomt.h5ad` | 03_3 | latent space, dimension stats (`vanished`), OOD/IND scores, UMAP |
| `shiao_tum_drvi_tum_64_nomt.h5ad` | 03_3 | the malignant object (all genes) with the latent coordinates |
| `shiao_tum.h5ad` | 03_2 | the malignant object, all genes, log-normalised: 03_6 scores on it |
| `shiao_tum_hvg_2k_nomt.h5ad` | 03_2 | the DRVI input, the ORA background of 03_7 |
| `cnv_status.csv` | 03_1 | the per-cell malignant call |
| `cell_annotation.csv` | 03_1 | CellTypist on the non-malignant cells |
| `drvi_leiden_<res>` columns of `embed_drvi_tum_64_nomt.h5ad` | 03_3 `run_drvi_tum.py` | the Leiden partitions of the DRVI states |

The last one matters as much as the model: the `STATE_PICKS` clusters are numbered on that
partition, and the kNN graph behind Leiden is not bit-reproducible across machines.

## Outputs

Heavy objects go to `$DATA_DIR/03_tum/` (see `datasets/README.md`). The repo keeps
`tables/<step>/` and `figures/<step>/`, one folder per step (03_8 adds `dr_umap/`, `state_de_*/` and
`state_dr_enrichment_*/`, 03_8b one folder per max(DR) variant); collection and run id are in each file name. Every table starts with a one-line header naming the run and
the pruning.
