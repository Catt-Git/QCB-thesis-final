# QCB-thesis-final

Code for my Master's thesis in Quantitative and Computational Biology (QCB): a single-cell
RNA-seq analysis of triple-negative breast cancer (TNBC) treated with pembrolizumab and
radiotherapy. The repository goes from the raw FASTQ files on SRA to a benchmark of
batch-correction methods, and then to an interpretable latent space (DRVI) of the malignant cells,
read against gene signatures of stemness, immunogenicity and EMT and followed across the
treatment timepoints.

> [!NOTE]
> The thesis PDF and the presentation slides (PDF and PPT) will be added here when they are ready.

## The dataset

Shiao et al., *Cancer Cell* 2024 (DOI 10.1016/j.ccell.2023.12.012), BioProject PRJNA1032700 /
GEO GSE246613. Each TNBC patient was biopsied three times: before treatment (`BASE`), after one
cycle of pembrolizumab (`PD1`) and after pembrolizumab plus radiotherapy (`RTPD1`). Patients are
grouped by response: `R1` early responders, `R2` responders only after radiotherapy, `NR`
non-responders. After processing: 149 libraries (CD45+ and CD45− sorts), 34 patients,
619,693 cells.

## The pipeline

Four phases, run in order; each folder has its own README with requirements, run order, design
choices and outputs.

| Phase | Folder | What it does | Main output |
|---|---|---|---|
| 00 | [00_raw_data_processing](00_raw_data_processing/README.md) | download, sample mapping, Cell Ranger | one `.h5ad` of raw counts with the sample metadata |
| 01 | [01_pre_processing](01_pre_processing/README.md) | QC, Scrublet, scran normalization, cell cycle, CellTypist annotation | `shiao.h5ad`: 619,693 cells × 30,869 genes, 48 cell types |
| 02 | [02_integration_benchmark](02_integration_benchmark/README.md) | ten batch-correction methods (batch = patient) scored with the scIB metrics | one summary table and the integrated embeddings |
| 03 | [03_drvi_tumoral_epi](03_drvi_tumoral_epi/README.md) | inferCNV malignant call, DRVI on the malignant cells, two routes of signature reading and their convergence, paired treatment tests | the named programs of the tumor cells and how they move with treatment |

Supporting folders:

- [environments](environments/README.md): the conda environments. `benchmark-py-r` runs almost
  everything; scGen, CytoTRACE2 and inferCNV have isolated environments.
- [datasets](datasets/README.md): what each phase reads and writes under `DATA_DIR`, plus the
  small files versioned in git (gene signatures, cell-cycle genes).

## Getting started

All heavy data (FASTQ, Cell Ranger outputs, `.h5ad` objects) live outside git, under one data
root that every script reads from the `DATA_DIR` environment variable:

```bash
conda env create -f environments/benchmark-py-r.yml
conda activate benchmark-py-r
export DATA_DIR=<data root>
```

Replace `<data root>` with your own path; to use the `datasets/` folder of this clone,
`export DATA_DIR="$(git rev-parse --show-toplevel)/datasets"`. The scripts stop with an error if
`DATA_DIR` is not set. Then follow the README of each phase, from 00 to 03.

Phase 00 needs a SLURM cluster and ~4.8 TB of disk (Cell Ranger on 149 libraries, ~2.4 TB of
FASTQ); phases 02 and 03 use a GPU for the
deep-learning models. Phase 03 can start from the checkpoints listed in its README instead of
retraining DRVI.

## Contact

For questions, or to get the already processed datasets instead of downloading and processing
everything again, write to Alberto Catalano, alberto.catalano02@outlook.it.

## Acknowledgements

I thank the laboratories and the PIs who hosted me during this thesis, for their guidance and
support.
