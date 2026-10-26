#!/usr/bin/env bash
#
# 03_1 infercnv_all: prepare the inputs and run inferCNV on every cohort.
#
#   1. prepare_infercnv_input.py   shiao.h5ad -> input/<cohort>/          [benchmark-py-r]
#   2. run_infercnv.R per cohort   input -> summary/<cohort>_cnv.csv      [infercnv-r]
#
# Cohorts run one at a time: inferCNV already uses --threads and holds the whole cohort in
# memory. Cohorts with a summary are skipped, so a re-run resumes; a failed cohort is
# reported at the end without stopping the others. The malignant call is made afterwards,
# in call_malignant.ipynb. --slurm submits the whole driver as one job (utils/slurm.sh),
# with one CPU per inferCNV thread.
#
# Usage:
#   export DATA_DIR=<data root>
#   bash infercnv_all.sh                              # prepare + every cohort, resuming
#   bash infercnv_all.sh --slurm --threads 16         # same, as one job on the cluster
#   bash infercnv_all.sh --force                      # re-run everything, overwriting
#   bash infercnv_all.sh --dry-run                    # print what would run, do nothing
#   bash infercnv_all.sh --threads 16                 # threads for infercnv::run() [default: 8]
#   bash infercnv_all.sh --leiden-resolution 0.01     # subcluster resolution [default: 0.005]
#   bash infercnv_all.sh --cohorts Patient52 Patient16  # only these cohorts
#   bash infercnv_all.sh prepare                      # only stage 1
#   bash infercnv_all.sh infercnv                     # only stage 2
#
# Environments can be overridden with PREP_ENV / INFERCNV_ENV.
# Logs go to logs/infercnv_all_<timestamp>.log as well as to the terminal (and to
# logs/infercnv_all_<jobid>.out with --slurm).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$SCRIPT_DIR/logs"

PREP_ENV="${PREP_ENV:-benchmark-py-r}"
INFERCNV_ENV="${INFERCNV_ENV:-infercnv-r}"
THREADS="${THREADS:-8}"
# Same verdict over 0.005-0.05 (sensitivity_resolution.sh); at the coarse end most cells
# sit in subclusters large enough to take a median over.
RESOLUTION="${RESOLUTION:-0.005}"

ARGS=("$@")
FORCE=0; DRY_RUN=0; HMM=0; USE_SLURM=0
COHORTS=()
STAGE_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --hmm) HMM=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --leiden-resolution) RESOLUTION="$2"; shift 2 ;;
    --cohorts) shift; while [ $# -gt 0 ] && [[ "$1" != -* ]]; do COHORTS+=("$1"); shift; done ;;
    -h|--help) sed -n '2,/^$/p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    prepare|infercnv) STAGE_FILTER+=("$1"); shift ;;
    *) echo "unknown stage '$1'; valid stages: prepare infercnv" >&2; exit 1 ;;
  esac
done

: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
CNV_DIR="$DATA_DIR/03_tum"
export FIG_DIR="${FIG_DIR:-$PHASE_DIR/figures/03_1_infercnv}"

if [ "$USE_SLURM" -eq 1 ]; then
  . "$PHASE_DIR/utils/slurm.sh"
  submit_driver infercnv_all "$THREADS" 64G "$LOG_DIR" "$DRY_RUN" "$SCRIPT_DIR/infercnv_all.sh" \
    ${ARGS[@]+"${ARGS[@]}"}
  exit 0
fi

run_stage() { [ ${#STAGE_FILTER[@]} -eq 0 ] || printf '%s\n' "${STAGE_FILTER[@]}" | grep -qx "$1"; }

# conda run, so the two environments do not leak into each other across the loop.
conda_run() { local env="$1"; shift; conda run --no-capture-output -n "$env" "$@"; }

mkdir -p "$CNV_DIR" "$LOG_DIR" "$FIG_DIR"
LOG_FILE="$LOG_DIR/infercnv_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "DATA_DIR    : $DATA_DIR"
echo "output dir  : $CNV_DIR"
echo "FIG_DIR     : $FIG_DIR"
echo "environments: $PREP_ENV (prepare) | $INFERCNV_ENV (inferCNV)"
echo "threads     : $THREADS | HMM: $HMM | leiden resolution: $RESOLUTION"
[ "$DRY_RUN" -eq 0 ] && echo "log         : $LOG_FILE"
echo
df -h "$DATA_DIR" | tail -1
echo

overall_start=$SECONDS

# Stage 1: prepare
if run_stage prepare; then
  echo "==================================================================="
  echo ">>> prepare  (prepare_infercnv_input.py)  $(date '+%Y-%m-%d %H:%M:%S')"
  echo "==================================================================="
  prep_args=()
  [ "$FORCE" -eq 1 ] && prep_args+=(--force)
  [ ${#COHORTS[@]} -gt 0 ] && prep_args+=(--cohorts "${COHORTS[@]}")
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] conda run -n $PREP_ENV python3 prepare_infercnv_input.py ${prep_args[*]:-}"
  else
    # ${a[@]+"${a[@]}"}: "${a[@]:-}" would pass an empty argument that argparse rejects.
    conda_run "$PREP_ENV" python3 "$SCRIPT_DIR/prepare_infercnv_input.py" ${prep_args[@]+"${prep_args[@]}"}
  fi
  echo
fi

# Stage 2: inferCNV per cohort
if run_stage infercnv; then
  CENSUS="$CNV_DIR/cohort_census.csv"
  [ -f "$CENSUS" ] || { echo "missing $CENSUS; run the prepare stage first" >&2; exit 1; }

  if [ ${#COHORTS[@]} -gt 0 ]; then
    TO_RUN=("${COHORTS[@]}")
  else
    # Column 1 = cohort, column 2 = status; only 'prepared' rows have inputs on disk.
    mapfile -t TO_RUN < <(awk -F, 'NR>1 && $2=="prepared" {print $1}' "$CENSUS")
  fi
  echo "cohorts to run: ${#TO_RUN[@]}"
  echo

  n_run=0; n_have=0; n_fail=0
  FAILED=()
  for cohort in "${TO_RUN[@]}"; do
    summary="$CNV_DIR/summary/${cohort}_cnv.csv"
    if [ "$FORCE" -eq 0 ] && [ -f "$summary" ]; then
      echo "[have] $cohort: $(basename "$summary") already exists, skipping"
      n_have=$((n_have + 1)); continue
    fi
    r_args=(--cohort "$cohort" --threads "$THREADS" --leiden-resolution "$RESOLUTION")
    [ "$FORCE" -eq 1 ] && r_args+=(--force)
    [ "$HMM" -eq 1 ] && r_args+=(--hmm)
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "[dry] conda run -n $INFERCNV_ENV Rscript run_infercnv.R ${r_args[*]}"
      n_run=$((n_run + 1)); continue
    fi

    echo "==================================================================="
    echo ">>> $cohort  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "==================================================================="
    step_start=$SECONDS
    # A failed cohort (e.g. degenerate reference) must not stop the others.
    if conda_run "$INFERCNV_ENV" Rscript "$SCRIPT_DIR/run_infercnv.R" "${r_args[@]}"; then
      echo "[ok] $cohort in $(( (SECONDS - step_start) / 60 ))m$(( (SECONDS - step_start) % 60 ))s"
      n_run=$((n_run + 1))
    else
      echo "[FAIL] $cohort (exit $?), continuing" >&2
      n_fail=$((n_fail + 1))
      FAILED+=("$cohort")
    fi
    # run_infercnv.R deletes its working directory only on success.
    [ -d "$CNV_DIR/work/$cohort" ] && rm -rf "$CNV_DIR/work/$cohort"
    df -h "$DATA_DIR" | tail -1
  done

  echo "==================================================================="
  echo "inferCNV: $n_run run, $n_have already done, $n_fail failed"
  [ "${#FAILED[@]}" -gt 0 ] && echo "failed cohorts: ${FAILED[*]}"
fi

echo "done in $(( (SECONDS - overall_start) / 60 ))m$(( (SECONDS - overall_start) % 60 ))s"
echo "next: call_malignant.ipynb on $CNV_DIR/summary/"
exit 0
