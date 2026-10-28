# 03_9_treatment — the named programs across the three timepoints

What did the treatment do to the tumor cells? Each patient was biopsied before treatment
(`BASE`), after one cycle of pembrolizumab (`PD1`) and after pembrolizumab plus radiotherapy
(`RTPD1`); patients are grouped by response (`NR`, `R1` early responders, `R2` responders only
after RT).

`treatment_programs_tum.py` tests whether a named program (a bar of `programs_named`) goes
up or down, paired within the patient. `treatment_common.py` holds the samples, the readouts and
the hypotheses.

```bash
export DATA_DIR=<data root>
python3 treatment_programs_tum.py
python3 treatment_programs_tum.py --min-cells 30
python3 treatment_programs_tum.py --no-gene-level   # skip the check that needs 03_6's scores
```

## One joint DRVI run

`drvi_tum_64_nomt` corrects the patient (`batch_key='cohort'`); the timepoint was never given to
the model, so the treatment effect stays in one shared coordinate system. One run per timepoint
would give axes with no correspondence, trained on different numbers of cells. Do not retrain
with `treatment` as a covariate: it would remove the signal measured here.

## The unit is the biopsy

Testing 42,000 cells would be pseudoreplication. The sample is `cohort × treatment` (CD45+ and
CD45− libraries pooled), with at least `--min-cells` malignant cells (default 20), and every
contrast is paired within the patient.

## The test

- **Readouts.** Each named direction (e.g. `DR 18-`) is oriented so higher means more of the
  program and z-scored over all cells, so deltas are in sd. A program readout is the mean of
  its directions.
- **Test.** Wilcoxon signed-rank on the per-patient difference of biopsy means, for
  `PD1_vs_BASE`, `RTPD1_vs_BASE` and `RTPD1_vs_PD1`. Fewer than 6 pairs: not tested.
- **Two BH families.** The hypotheses fixed before looking at the data (`HYPOTHESES`):

  | program | contrast | expected |
  |---|---|---|
  | HALLMARK_IFNG, MP17_INTERFERON_MHC_II_1 | PD1 vs BASE | up (IFN-γ from re-activated T cells) |
  | HALLMARK_IFNA | RTPD1 vs PD1 | up (RT-induced type-I IFN) |
  | MP19_EPITHELIAL_SENESCENCE | RTPD1 vs PD1 | up (RT-induced senescence) |
  | MP1 / MP2 cell cycle | RTPD1 vs PD1 | down (RT-induced arrest) |

  Everything else is one exploratory family; a hit there is a lead.
- **Gene-level check** (a priori programs only): the same test on 03_6's raw signature score.
  A DR that moves while its genes do not is a property of the embedding.
- **Response** is described, not tested as an interaction (2–8 patients per group and
  contrast): median delta per group with its n, plus an exploratory Kruskal–Wallis.

## Inputs

| file | written by | where |
|---|---|---|
| `embed_drvi_tum_64_nomt.h5ad` | 03_3 | `$DATA_DIR/03_tum/` |
| `program_dimensions_consensus_drvi_tum_64_nomt.csv` | 03_8 `dr_naming_tum.py` | `tables/03_8_convergence/` |
| `signature_scores_{scie,gavish_tnbc}_drvi_tum_64_nomt.csv` (optional, gene-level check) | 03_6 | `$DATA_DIR/03_tum/` |

## Outputs

`tables/03_9_treatment/`:

- `treatment_cell_counts_*`: malignant cells per patient × timepoint
- `treatment_pseudobulk_*`: biopsy means and the fraction above 1 sd, per readout
- `treatment_paired_tests_*`: one row per readout × contrast (family, FDR, per-response medians, gene-level columns)

`figures/03_9_treatment/`:

- `a_priori_trajectories`: one line per patient, response medians in bold
- `program_trajectories`: the same for every named program
- `paired_delta_programs`, `paired_delta_directions`: median paired change, a priori cells
  outlined; `*` FDR < 0.05, `°` p < 0.05 before FDR (a lead)

## On `drvi_tum_64_nomt`

20 patients, 49 biopsies with ≥ 20 malignant cells: 18 pairs for `PD1_vs_BASE`, 11 for
`RTPD1_vs_BASE`, 10 for `RTPD1_vs_PD1`. 22 programs and their 30 directions, 156 tests (6 a
priori, 150 exploratory).

| program | contrast | median Δ (sd) | patients | FDR | gene level p |
|---|---|---|---|---|---|
| MP1_CELL_CYCLE_G2_M | RTPD1 vs PD1 | −0.36 | 8 of 10 down | 0.029 | 0.004 (FDR 0.023) |
| MP2_CELL_CYCLE_G1_S | RTPD1 vs PD1 | −0.30 | 8 of 10 down | 0.029 | 0.19 |
| MP19_EPITHELIAL_SENESCENCE | RTPD1 vs PD1 | +0.41 | 8 of 10 up | 0.13 | 0.037 |
| HALLMARK_IFNA | RTPD1 vs PD1 | +0.18 | 8 of 10 up | 0.13 | 0.049 |
| HALLMARK_IFNG | PD1 vs BASE | +0.11 | 12 of 18 up | 0.13 | 0.043 |
| MP17_INTERFERON_MHC_II_1 | PD1 vs BASE | +0.08 | 12 of 18 up | 0.21 | 0.043 |

All six move the expected way; only the radiotherapy arrest of the cell cycle passes, and only
G2/M holds on the genes themselves (the other gene-level FDRs are 0.059). No exploratory test
reaches FDR < 0.05 (lowest 0.098); 39 of 150 have p < 0.05, 26 of them on `RTPD1_vs_BASE`.

## Limits

- **One biopsy per timepoint**: the timepoint cannot be separated from that biopsy's library
  or sampling site.
- **Survivor bias**: in responders the late malignant cells are the survivors, so a shift can be
  selection as much as reprogramming.
- **Power**: 20 patients, 10–18 pairs per contrast, 2–8 per response group; a null result is
  weak evidence.
- **Programs named on all timepoints together**: a state that appears only after treatment
  and is too rare to earn a direction is invisible here.
