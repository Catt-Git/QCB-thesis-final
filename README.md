# QCB-thesis-final
Final repository for my Master's thesis, containing the thesis PDF, presentation slides  and curated, organized code used for the analyses and results presented in the thesis.


Contattatemi per avere i dataset già processati in modo da non dover riscaricare tutto, o per qualsiasi altra domanda

## Before running anything: DATA_DIR

All heavy data (FASTQ, Cell Ranger outputs, .h5ad objects) live outside git, under a single data root that every script reads from the DATA_DIR environment variable. The Usage sections of the scripts and READMEs write it as a placeholder:

```bash
export DATA_DIR=<data root>
```

Replace `<data root>` with your own path before running any command. If you keep the data in the datasets/ folder of this clone, you can use:

```bash
export DATA_DIR="$(git rev-parse --show-toplevel)/datasets"
```

The scripts stop with an error if DATA_DIR is not set. See [datasets/README.md](datasets/README.md) for what each phase writes there.
