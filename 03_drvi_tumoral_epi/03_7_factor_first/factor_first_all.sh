#!/usr/bin/env bash
#
# 03_7 factor-first (Route B): the enrichment on every collection, in order:
#
#   1. factor_first_tum.py --collection scie    -> dim_geneset_signed_significance_scie + figure
#   2. factor_first_tum.py --collection emt     -> same, on the EMT lists
#   3. factor_first_tum.py --collection gavish  -> same, on the TNBC metaprograms
#
# Needs 03_4 (the .gmt), 03_3 (the embedding) and 03_6 (each collection's row order): run
# 03_6_cell_first/cell_first_all.sh first. Light (no model, no GPU, no network), so it runs
# locally. Steps with an existing output are skipped, so a re-run resumes; --force re-runs
# everything (--overwrite). --slurm submits it as one job.
#
# Usage:
#   export DATA_DIR=<data root>
#   ./factor_first_all.sh                     # every collection, resuming
#   ./factor_first_all.sh --slurm             # same chain as one job on the cluster
#   ./factor_first_all.sh --force             # re-run everything, overwriting
#   ./factor_first_all.sh --dry-run           # print what would run, do nothing
#   ./factor_first_all.sh --all-metaprograms  # gavish on all 41 metaprograms
#   ./factor_first_all.sh emt                 # only the named step(s), in file order
#
# Step names: scie, emt, gavish.
# Logs go to 03_7_factor_first/logs/factor_first_all_<timestamp>.log and to the terminal (and
# to factor_first_all_<jobid>.out with --slurm).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"

STEP_NAMES=(scie emt gavish)

ARGS=("$@")
FORCE=0; DRY_RUN=0; USE_SLURM=0; ALL_METAPROGRAMS=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
    --all-metaprograms) ALL_METAPROGRAMS=1; shift ;;
    -h|--help) sed -n '3,/^$/p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

for requested in "${STEP_FILTER[@]:-}"; do
  [ -n "$requested" ] || continue
  in_list "$requested" "${STEP_NAMES[@]}" || {
    echo "unknown step '$requested'; valid steps: ${STEP_NAMES[*]}" >&2; exit 1; }
done

# Before resolve_run.py: the submitting shell need not have the Python env.
if [ "$USE_SLURM" -eq 1 ]; then
  . "$PHASE_DIR/utils/slurm.sh"
  submit_driver factor_first_all 2 16G "$LOG_DIR" "$DRY_RUN" "$SCRIPT_DIR/factor_first_all.sh" \
    ${ARGS[@]+"${ARGS[@]}"}
  exit 0
fi

# Run id and output root from Python, so the resume check looks where the scripts write.
read -r RUN_ID COMPARTMENT OUT_ROOT_DIR < <(cd "$PHASE_DIR" && "$PYTHON" utils/resolve_run.py)
TABLE_ROOT="$OUT_ROOT_DIR/tables"
if [ "$OUT_ROOT_DIR" != "$PHASE_DIR" ]; then LOG_DIR="$OUT_ROOT_DIR/logs"; fi

# Mirrors sig_collections.resolve(): gavish writes under gavish_tnbc unless --all-metaprograms.
GAVISH_SLUG="gavish_tnbc"
[ "$ALL_METAPROGRAMS" -eq 1 ] && GAVISH_SLUG="gavish"
SLUGS=(scie emt "$GAVISH_SLUG")

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/factor_first_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "03_7 factor-first"
echo "run id    : $RUN_ID  (compartment $COMPARTMENT)"
echo "DATA_DIR  : $DATA_DIR"
echo "gavish    : $GAVISH_SLUG"
[ "$DRY_RUN" -eq 0 ] && echo "log       : $LOG_FILE"
echo

n_run=0; n_have=0
for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  slug="${SLUGS[$i]}"
  output="$TABLE_ROOT/03_7_factor_first/dim_geneset_signed_significance_${slug}_${RUN_ID}.csv"
  row_order="$TABLE_ROOT/03_6_cell_first/dimension_row_order_${slug}_${RUN_ID}.csv"

  if [ ${#STEP_FILTER[@]} -gt 0 ] && ! in_list "$name" "${STEP_FILTER[@]}"; then
    continue
  fi

  if [ "$FORCE" -eq 0 ] && [ -e "$output" ]; then
    echo "[have] $name  -> $output"
    n_have=$((n_have + 1)); continue
  fi

  cmd=("$PYTHON" factor_first_tum.py --collection "$name")
  # An `if`, not `cond && ...`, which would trip set -e.
  if [ "$ALL_METAPROGRAMS" -eq 1 ] && [ "$name" = "gavish" ]; then cmd+=(--all-metaprograms); fi
  if [ "$FORCE" -eq 1 ]; then cmd+=(--overwrite); fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry ] $name  -> ${cmd[*]}"
    n_run=$((n_run + 1)); continue
  fi

  [ -e "$row_order" ] || {
    echo "[error] $name needs 03_6's $row_order; run 03_6_cell_first/cell_first_all.sh first" >&2
    exit 1; }

  echo "[run ] $name  -> ${cmd[*]}"
  ( cd "$SCRIPT_DIR" && "${cmd[@]}" )
  [ -e "$output" ] || { echo "[error] $name finished but $output is missing" >&2; exit 1; }
  echo "[done] $name"
  echo
  n_run=$((n_run + 1))
done

echo "done: $n_run step(s) run, $n_have already done"
