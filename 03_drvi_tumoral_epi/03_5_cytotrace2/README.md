# 03_5_cytotrace2

CytoTRACE2 potency per cell, one patient at a time. A Route A readout with no Route B
counterpart: independent of the gene lists, it is the only non-circular evidence on the
stemness axis.

```bash
conda activate cytotrace2-py          # not benchmark-py-r, see below
export DATA_DIR=<data root>
python3 cytotrace2_tum.py --dry-run   # export the per-patient matrices, do not score
python3 cytotrace2_tum.py             # export + score + concatenate -> cytotrace2_drvi_tum_64_nomt.csv
python3 cytotrace2_tum.py --cores 4   # cap the cores (1-2 under 16 GB)
```

The score is computed on raw counts and never sees the latent space; the run id in the file
name is only how 03_6 finds it.

## Its own environment

`cytotrace2-py` pins `numpy<2.0.0`: installed into `benchmark-py-r` it would roll back numpy,
pandas, scanpy and anndata. Use `environments/cytotrace2-py.yml`; the script checks for the
package and stops otherwise. Without its output `cell_first_tum.py` drops the CytoTRACE2
readout and says so.

## Why per patient

The score is relative to the object it is given, so scoring all cohorts together would let
it absorb batch (e.g. depth). Each cohort is scored on its own, so scores are comparable within
a patient and only ordinally across patients, which is how 03_6 uses them (z within `cohort`).

## Potency or proliferation?

Every cell is aneuploid, so a high score cannot be a normal progenitor; it can still be a
cycling cell (`Lumsec-prol` is 53% of the subset). The summary tables are therefore per
`drvi_leiden`, per pre-CNV CellTypist label and per `phase`: if the top of the Leiden table is
the S/G2M clusters and the top of the label table is `Lumsec-prol`, potency is tracking the
cycle.

## What comes out

```
$DATA_DIR/03_tum/cytotrace2_drvi_tum_64_nomt.csv          per cell: score, potency, relative, preKNN
$DATA_DIR/03_tum/cytotrace2_work_drvi_tum_64_nomt/        per-patient matrices and raw results
../tables/03_5_cytotrace2/cytotrace2_by_<grouping>_scie_drvi_tum_64_nomt.csv
```

The summary tables carry `scie`, the collection that uses CytoTRACE2. `--write-obs` also adds
the columns to `shiao_tum_drvi_tum_64_nomt.h5ad` (off by default: that object is 03_3's output).
The per-patient matrices are large and deleted afterwards unless `--keep-work`.
