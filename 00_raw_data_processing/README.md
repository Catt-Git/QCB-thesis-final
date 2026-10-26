# 00_raw_data_processing

First phase of the thesis: from the raw FASTQ files on SRA to a single `.h5ad` of raw counts,
one sample per library, carrying each sample's metadata (cohort, treatment, response, ...).

Dataset: Shiao et al., *Cancer Cell* 2024 (DOI 10.1016/j.ccell.2023.12.012),
BioProject **PRJNA1032700** / GEO **GSE246613**.

## Requirements

The whole phase runs on a SLURM cluster: the ~2.4 TB of FASTQ are downloaded straight onto its
storage and never leave it.

- Conda environment `benchmark-py-r` (`environments/benchmark-py-r.yml`, which also provides the SRA Toolkit).
- Cell Ranger 4.0.0 and its reference `refdata-cellranger-GRCh38-3.0.0` (~11 GB, the same as the
  paper), downloaded once and pointed to by `REF`:
  ```bash
  curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
  tar -xzf refdata-cellranger-GRCh38-3.0.0.tar.gz # REF = the extracted refdata-cellranger-GRCh38-3.0.0/
  ```
- ~4.8 TB free under `DATA_DIR` (the FASTQ plus their renamed copies).
- Outbound internet from the compute nodes (steps 1 and 3).

The `.slurm` files are written for the thesis cluster: their `module load` lines are
site-specific. `DATA_DIR` and `REF` have no default and must be exported before launching.

## Quick start

From the cluster's login node:

```bash
cd 00_raw_data_processing
export DATA_DIR=<data root>
export REF=<path to refdata-cellranger-GRCh38-3.0.0>
bash rawprocessing_all.sh --dry-run # what is already done, what would be submitted
bash rawprocessing_all.sh # submit the missing steps as one afterok chain
```

`rawprocessing_all.sh` submits the jobs below, each starting only if the previous one succeeded,
and skips what is already on disk, so running it again after a failure resumes where it stopped.
A single step can also be submitted by hand, from its own folder
(e.g. `cd 00_3_download && sbatch submit_download.slurm`).

## Execution order

| # | Folder | Script | What it does | Partition, CPUs, memory |
|---|--------|--------|--------------|-------------------------|
| 1 | `00_1_metadata` | `fetch_metadata.sh` | Download the ENA metadata of the BioProject | `normal`, 1, 2G |
| 2 | `00_2_sample_mapping` | `build_sample_map.sh` | Map SRR to patient/timepoint/fraction; drop TCR/BCR and mouse libraries, and runs with no FASTQ on ENA | `normal`, 1, 2G |
| 3 | `00_3_download` | `download_fastq.sh` | Download the FASTQ (`prefetch` + `fasterq-dump`), one array task per run | `long`, 8, 16G per task |
| 3b | `00_3_download` | `verify_integrity.sh` | Check every FASTQ: `gzip -t`, R1/R2 read lengths, 5' TSO | `long`, 2, 8G |
| 4 | `00_4_rename` | `rename_for_cellranger.sh` | Copy to `<sample>_S1_L001_R[12]_001.fastq.gz`; write `unique_samples.txt` | `long`, 2, 8G |
| 5 | `00_5_cellranger` | `run_cellranger_batch.sh` | `cellranger count`, 5 samples in parallel | `long`, 90, 470G |
| 6 | `00_6_build_h5ad` | `build_requested_metadata.py` + `build_combined_h5ad.py` | Build the per-sample metadata table, then merge the matrices into one AnnData | `normal`, 8, 150G |
| - | `utils` | `check_download_completeness.sh` | List the runs missing R1 or R2 | login node |
| - | `utils` | `check_cellranger_completeness.sh` | List the samples missing their Cell Ranger matrix | login node |

Each step folder holds the script that does the work and the `submit_*.slurm` that launches it.

## Data layout

Everything large lives outside the repo, under `$DATA_DIR`:

    $DATA_DIR/
    ├── fastq_raw/                 # step 3: SRR*_1/_2.fastq.gz as downloaded (~2.4 TB)
    ├── renamed_fastq/             # step 4: <sample>_S1_L001_R[12]_001.fastq.gz + unique_samples.txt
    ├── cellranger_out/            # step 5: <sample>/outs/filtered_feature_bc_matrix.h5
    └── all_samples_combined.h5ad  # step 6: final AnnData

The small tables (`ena_full_metadata.tsv`, `sample_map_gex_final.tsv`, `sample_metadata_final.tsv`)
are versioned in the repo, next to the script that writes them.

**Freeing space.** Once `all_samples_combined.h5ad` exists, `fastq_raw/` and `renamed_fastq/`
(~4.8 TB together, plus any `.sra` file `prefetch` left in `fastq_raw/`) are no longer read by
anything and can be deleted. If space is not a problem, keep everything until the analyses are
finished: rebuilding the FASTQ means downloading 2.4 TB again. `rawprocessing_all.sh` treats the phase as complete
while the `.h5ad` exists, so deleting the FASTQ does not trigger a new download.

## Key parameters

- **Download:** `prefetch --max-size 100G` (the 20 GB default skips the largest runs), then
  `fasterq-dump --include-technical --split-files`. Without `--include-technical` the barcode read
  (R1) is dropped as a technical read.
- **Cell Ranger:** version 4.0.0, reference GRCh38-3.0.0 (as in the paper),
  `--chemistry=SC5P-R2` (10x Chromium 5' v1, confirmed by the TSO `AAGCAGTGGTATCAAC` at the start
  of R2), `--r1-length=26` (16 bp barcode + 10 bp UMI; R1 is sequenced at 101 bp, 151 bp for two runs).
- **AnnData:** built from `filtered_feature_bc_matrix.h5` (Cell Ranger's cell calling). Raw counts
  in `.X` and in `.layers['counts']`; gene symbols made unique as `var_names`, Ensembl IDs in
  `.var['gene_ids']`; barcodes prefixed with the sample name; `ad.concat(join='outer', merge='same')`.

## Expected output

`all_samples_combined.h5ad` should contain:

- **149 samples** (103 CD45+, 46 CD45-): the BioProject has 150 human GEX runs, and step 2 drops
  `SRR26541168` (`P30_C_P`), which has no FASTQ on ENA.
- **34 patients**, with response 9 R1 / 14 R2 / 11 NR.
- **`.obs` columns:** `sample`, `batch`, `cohort`, `treatment`, `fraction`, `dataset_origin`, `response`.
  - `fraction` and `dataset_origin` are per library: CD45+ (immune) or CD45- (non-immune), from
    the magnetic-bead sort.
  - `response` is per patient, transcribed from GSE246613 and not recomputed.

## Repairing a single run

If step 3b reports a file as `CORRUPTED`, delete that run's FASTQ and launch the wrapper again: it
downloads only what is missing and re-runs the steps after it.

```bash
rm -f "$DATA_DIR"/fastq_raw/SRR26540980_[12].fastq.gz
rm -rf "$DATA_DIR"/cellranger_out/P07_A_P # only if the run already went through Cell Ranger
rm -f "$DATA_DIR"/all_samples_combined.h5ad # only if it was already built
bash rawprocessing_all.sh
```

Each deletion is needed: the download skips a run whose FASTQ exist, Cell Ranger skips a sample
whose filtered matrix exists, and the wrapper considers the phase complete while the `.h5ad` exists. In the thesis run this happened to `P07_A_P` (`SRR26540980`)
and `P10_C_P` (`SRR26541023`). Their new copies have 151 bp R1 reads, which `--r1-length=26` handles.
