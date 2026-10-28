# 03_4_signatures

The prior knowledge before either route touches it. Reads the collaborator's plain-text lists
from `$DATA_DIR/signatures/` and writes one `.gmt` per collection, with the provenance in the
description field, so the same file feeds Route A and Route B.

## Requirements

- Conda environment `benchmark-py-r`.
- The gene lists in `$DATA_DIR/signatures/` (versioned in `datasets/signatures/`).
- `$DATA_DIR/03_tum/shiao_tum.h5ad` and `shiao_tum_hvg_2k_nomt.h5ad` from 03_2: the gene universe
  of Route A and the HVG background of Route B. Both are read backed, so the step is light.

## Quick start

```bash
cd 03_drvi_tumoral_epi/03_4_signatures
export DATA_DIR=<data root>
bash signatures_all.sh                               # all three collections + the EMT check, resuming
bash signatures_all.sh --slurm                       # the same as one job on the cluster
bash signatures_all.sh --force                       # re-run everything
bash signatures_all.sh --dry-run                     # print what would run
```

Existing outputs are skipped, so after a code change re-run with `--force`.
`signature_composition_tum.py` is not in the wrapper: it needs the 03_6 scores, so
`../03_6_cell_first/cell_first_all.sh` runs it after 03_6.

## The files

| File | What it does |
|---|---|
| `signatures_all.sh` | runs `build_signatures_tum.py` on scie, emt and gavish, then `emt_vs_gavish_tum.py`; `--slurm` for the cluster |
| `build_signatures_tum.py` | lists → `.gmt`, coverage, Jaccard, the two figures |
| `emt_vs_gavish_tum.py` | do the collaborator's EMT lists and Gavish's EMT metaprograms name the same genes? |
| `signature_composition_tum.py` | which genes carry each score, and are they measurable. Runs after 03_6, to join the real scores |

One script at a time:

```bash
python3 build_signatures_tum.py                      # scie, the default
python3 build_signatures_tum.py --collection emt     # the EMT lists
python3 build_signatures_tum.py --collection gavish  # the 22 TNBC-relevant metaprograms
python3 build_signatures_tum.py --collection gavish --all-metaprograms  # all 41 instead
python3 build_signatures_tum.py --allow-low-coverage # report, do not stop

python3 emt_vs_gavish_tum.py                     # the EMT lists against MP12 - MP15
python3 emt_vs_gavish_tum.py --all-metaprograms  # MP16 (glioma mesenchymal) too
```

By default the script stops when a list maps under 60% of its symbols to the object;
`--allow-low-coverage` prints a warning instead and writes the `.gmt` and the tables anyway.

## The gene universes

- Route A scores on all genes of the object: 25,197 on 41,943 malignant cells.
- Route B only tests the part of a signature inside the 2,000 HVGs DRVI was trained on: the
  `n_in_hvg_background` column of the coverage table says how much of each list that is.

## What comes out

`../tables/03_4_signatures/`:

```
signatures_<collection>_tum.gmt                  the collection as used (mapped genes)
coverage_<collection>_drvi_tum_64_nomt.csv       per list: mapped fraction, HVG count
jaccard_<collection>_drvi_tum_64_nomt.csv        pairwise overlap
shared_genes_<collection>_drvi_tum_64_nomt.csv   the same as counts
signature_concentration_<collection>_*.csv       signature_composition_tum.py: how many genes each score behaves as if it had
signature_gene_contribution_<collection>_*.csv   the same, per gene
```

Figures in `../figures/03_4_signatures/`. The `.gmt` is named after the object, not the run: it
depends only on the genes of `shiao_tum.h5ad`.

- **coverage**: an unmapped symbol is not measured, which is not the same as not expressed.
  Below `MIN_MAPPED_FRACTION = 0.60` the script stops (`--allow-low-coverage` to override).
- **jaccard**: the immune lists are largely nested and the embryonic stemness lists overlap
  heavily, so they are not independent tests. Not written for `emt`: one list per axis, nothing
  to compare.

`EMT_B_EPITHELIAL` (4 genes) and `EMT_B_HYBRID` (6) have fewer than 10 genes in the HVG
background, so they are weak on Route B; Route A still scores them in full.

`gavish_tnbc` maps well (lowest `MP7_STRESS_IN_VITRO`, 0.880; `MP1_CELL_CYCLE_G2_M` 50/50), but
five metaprograms have fewer than 10 genes in the HVG panel:

| metaprogram | genes in HVG background |
|---|---|
| `MP4_CHROMATIN` | 3 |
| `MP8_PROTEASOMAL_DEGRADATION` | 6 |
| `MP9_UNFOLDED_PROTEIN_RESPONSE` | 3 |
| `MP20_MYC` | 5 |
| `MP21_RESPIRATION` | 3 |

A low count limits power, it does not make the list untestable: `MP21_RESPIRATION` is still
convergent on `DR 3+` (ρ 0.389, FDR 8.7e-03), while `MP4_CHROMATIN` is the one metaprogram
Route B never finds. These are housekeeping-like programs, expressed everywhere and therefore
not highly variable.

## EMT twice: the collaborator's lists against Gavish's metaprograms

`emt_vs_gavish_tum.py`. The `emt` collection (collaborator list B: epithelial, hybrid,
mesenchymal; a hypothesis) and the EMT axis of `gavish_tnbc` (MP12 - MP15, a vocabulary) meet
in 03_8; this checks they are not the same genes twice.

```
jaccard_emt_vs_gavish_emt_drvi_tum_64_nomt.csv        lists × metaprograms
shared_genes_emt_vs_gavish_emt_drvi_tum_64_nomt.csv   the same as counts
jaccard_bands_emt_vs_gavish_emt_drvi_tum_64_nomt.csv  the reference bands
emt_vs_gavish_pairs_emt_drvi_tum_64_nomt.csv          the cross pairs, with the shared genes
```

The figure (`jaccard_emt_vs_gavish_emt_*.png`) has two panels: on the left the whole matrix, on
the right the lists × metaprograms block on its own scale, with the shared-gene counts.

They are not: the lists overlap the metaprograms about as little as the metaprograms overlap
each other. The small overlap there is is still structured, one metaprogram per axis:

| axis | best MP | shared genes |
|---|---|---|
| mesenchymal | MP12_EMT_1 | ACTA2, FN1, MMP2, SNAI2, SPARC, VIM |
| hybrid | MP13_EMT_2 | ITGAV, LAMC2, PDPN, TNC |
| epithelial | MP14_EMT_3 | CD24, CLDN4, KRT8, KRT18, KRT19, LMNA |
| — | MP15_EMT_4 | nothing above 0.03 |

MP15 matches no list: a dimension whose top metaprogram is MP15 cannot be checked against the
collaborator's lists.
