#!/usr/bin/env bash
#
# 03_2 subsetting: the scripts that turn the raw subset of subset_and_qc.ipynb
# into the definitive object, in order:
#
#   1. scran_norm_tum.py       shiao_tum_raw.h5ad     -> shiao_tum_norm.h5ad
#   2. cell_cycle_score_tum.py shiao_tum_norm.h5ad    -> shiao_tum_norm_cc.h5ad
#   3. reduce_data_tum.py      shiao_tum_norm_cc.h5ad -> shiao_tum.h5ad
#                                                     +  shiao_tum_hvg_2k_nomt.h5ad (+ _list.csv)
#
# ~42k cells, so the chain runs locally in sequence; --slurm submits it as one job
# (utils/slurm.sh). Steps with an existing output are skipped, so a re-run resumes;
# delete an .h5ad truncated by a crash first. --force re-runs everything.
#
# Usage:
#   export DATA_DIR=<data root>
#   ./subsetting_all.sh                 # every step, resuming
#   ./subsetting_all.sh --slurm         # same chain as one job on the cluster
#   ./subsetting_all.sh --force         # re-run everything, overwriting
#   ./subsetting_all.sh --dry-run       # print what would run, do nothing
#   ./subsetting_all.sh cc reduce       # only the named step(s), in file order
#
# Step names: norm, cc, reduce.
# Logs go to 03_2_subsetting/logs/subsetting_all_<timestamp>.log and to the terminal (and to
# subsetting_all_<jobid>.out with --slurm).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"

PREFIX="shiao_tum"

STEP_NAMES=(norm cc reduce)
STEP_SCRIPTS=(
  scran_norm_tum.py
  cell_cycle_score_tum.py
  reduce_data_tum.py
)
STEP_SUFFIXES=(
  _norm.h5ad
  _norm_cc.h5ad
  .h5ad
)

ARGS=("$@")
FORCE=0; DRY_RUN=0; USE_SLURM=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
    -h|--help) sed -n '3,/^$/p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
TUM_DIR="$DATA_DIR/03_tum"

# The Tirosh/Regev list needed by the cc step. Same default as 01_4.
export CC_GENES="${CC_GENES:-$DATA_DIR/regev_lab_cell_cycle_genes.txt}"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

for requested in "${STEP_FILTER[@]:-}"; do
  [ -n "$requested" ] || continue
  in_list "$requested" "${STEP_NAMES[@]}" || {
    echo "unknown step '$requested'; valid steps: ${STEP_NAMES[*]}" >&2; exit 1; }
done

# The input of the whole chain comes from the notebook, not from this script.
INPUT_H5AD="$TUM_DIR/${PREFIX}_raw.h5ad"
if [ ! -f "$INPUT_H5AD" ] && [ "$DRY_RUN" -eq 0 ]; then
  echo "missing $INPUT_H5AD" >&2
  echo "run subset_and_qc.ipynb first" >&2
  exit 1
fi

if [ "$USE_SLURM" -eq 1 ]; then
  . "$(dirname "$SCRIPT_DIR")/utils/slurm.sh"
  submit_driver subsetting_all 8 64G "$LOG_DIR" "$DRY_RUN" "$SCRIPT_DIR/subsetting_all.sh" \
    ${ARGS[@]+"${ARGS[@]}"}
  exit 0
fi

mkdir -p "$TUM_DIR" "$LOG_DIR"
LOG_FILE="$LOG_DIR/subsetting_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "DATA_DIR  : $DATA_DIR"
echo "prefix    : $PREFIX"
echo "output dir: $TUM_DIR"
echo "CC_GENES  : $CC_GENES"
echo "python    : $($PYTHON --version 2>&1)"
[ "$DRY_RUN" -eq 0 ] && echo "log       : $LOG_FILE"
echo

n_run=0; n_have=0; n_skip=0
overall_start=$SECONDS

for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  script="$SCRIPT_DIR/${STEP_SCRIPTS[$i]}"
  output="$TUM_DIR/${PREFIX}${STEP_SUFFIXES[$i]}"

  if [ ${#STEP_FILTER[@]} -gt 0 ] && ! in_list "$name" "${STEP_FILTER[@]}"; then
    n_skip=$((n_skip + 1)); continue
  fi

  if [ "$FORCE" -eq 0 ] && [ -f "$output" ]; then
    echo "[have] $name: $output already exists, skipping"
    n_have=$((n_have + 1)); continue
  fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] $name: $PYTHON $script -> $output"
    n_run=$((n_run + 1)); continue
  fi

  echo "==================================================================="
  echo ">>> $name  ($(basename "$script"))  $(date '+%Y-%m-%d %H:%M:%S')"
  echo "==================================================================="
  step_start=$SECONDS
  ( cd "$SCRIPT_DIR" && "$PYTHON" "$script" )

  # Stop before the next step reads a stale file.
  [ -f "$output" ] || { echo "[error] $name finished but $output is missing" >&2; exit 1; }

  echo "[ok] $name in $(( (SECONDS - step_start) / 60 ))m$(( (SECONDS - step_start) % 60 ))s -> $output"
  n_run=$((n_run + 1))
done

echo "==================================================================="
echo "done in $(( (SECONDS - overall_start) / 60 ))m$(( (SECONDS - overall_start) % 60 ))s: \
$n_run step(s) run, $n_have already done, $n_skip filtered out"
[ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
  echo "      (--force to re-run and overwrite the $n_have existing output(s))"
echo "next: visualization_tum.ipynb on $TUM_DIR/${PREFIX}.h5ad"
exit 0
