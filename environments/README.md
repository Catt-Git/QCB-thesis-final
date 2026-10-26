# environments

Conda/mamba environment files used across the thesis.

## The environments

| File | Env name | Where / when to use |
|---|---|---|
| `benchmark-py-r.yml` | `benchmark-py-r` | **Main environment.** Python 3.12 + GPU (PyTorch 2.13, CUDA 12.9) + the R/Bioconductor stack (scran, batchelor, rpy2, anndata2ri) |
| `scgen-py.yml` | `scgen-py` | **Isolated env for scGen only.** Old pinned deps (Python 3.9, `anndata 0.10`, `scanpy 1.9`, `torch 2.0`, `scvi-tools 0.19.0`, `scgen 2.1.0`) that conflict with the main stack, so scGen gets its own env. |
| `cytotrace2-py.yml` | `cytotrace2-py` | **Isolated env for CytoTRACE2 only** (`03_drvi_tumoral_epi/03_5_cytotrace2`). `cytotrace2-py` pins `numpy<2.0.0`, so it cannot live beside the main stack; CPU-only torch, Python 3.11. |
| `infercnv-r.yml` | `infercnv-r` | **Isolated env for inferCNV only** (`03_drvi_tumoral_epi/03_1_infercnv`). R 4.5 + `bioconductor-infercnv` + JAGS. R only (the Python half of that step runs in `benchmark-py-r`). |

## Usage

Create an environment from its file (use `mamba` if available, it's faster):

```bash
conda env create -f environments/benchmark-py-r.yml
conda activate benchmark-py-r
```

Then run the code as described in each phase's README (`export DATA_DIR=<data root>` first).

- Local machine: `benchmark-py-r`.
- Cluster / SLURM: `benchmark-py-r` (same env, same name).
- scGen integration method only: `scgen-py`.
- CytoTRACE2 (03_5) only: `cytotrace2-py`.
- inferCNV (03_1) only: `infercnv-r`. That step runs its Python half in `benchmark-py-r` and
  switches environment per stage; `03_1_infercnv/infercnv_all.sh` does it with `conda run`.

## R inside the metrics environment

No metric of the 12 in `02_integration_benchmark` is computed through R, but the
metrics code imports `anndata2ri`, which starts an embedded R: the environment
still has to be *activated* (not just its `bin/python` invoked directly) so that
rpy2 can find R. `02_integration_benchmark/utils/smoke_test_metrics.py` prints
the R version it reaches, which is the cheap way to check this before launching
anything long.
