#!/usr/bin/env bash
#
# 03_1 sensitivity_reference: prepare, run and analyze the four reference configurations
# on two cohorts (8 inferCNV runs, ~60 min). Writes only under
# $DATA_DIR/03_tum/sensitivity/reference/ and ../tables/03_1_infercnv/sensitivity/.
#
# Usage:
#   export DATA_DIR=<data root>
#   bash utils/sensitivity_reference.sh
#   THREADS=8 bash utils/sensitivity_reference.sh    # threads per run [20]
#   bash utils/sensitivity_reference.sh --slurm      # same, as one job on the cluster

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS="${THREADS:-20}"

case "${1:-}" in
  "") ;;
  --slurm)
    : "${DATA_DIR:?set DATA_DIR}"
    STEP_DIR="$(dirname "$SCRIPT_DIR")"
    . "$(dirname "$STEP_DIR")/utils/slurm.sh"
    submit_driver sensitivity_reference "$THREADS" 64G "$STEP_DIR/logs" 0 \
      "$SCRIPT_DIR/sensitivity_reference.sh"
    exit 0 ;;
  -h|--help) sed -n '2,/^$/p' "${BASH_SOURCE[0]}"; exit 0 ;;
  *) echo "unknown option: $1" >&2; exit 1 ;;
esac

: "${DATA_DIR:?set DATA_DIR}"
ROOT="$DATA_DIR/03_tum/sensitivity/reference"

conda run --no-capture-output -n "${PREP_ENV:-benchmark-py-r}" \
  python3 "$SCRIPT_DIR/sensitivity_reference.py" --mode prepare || exit 1

for cfg in "$ROOT"/*/*/; do
  [ -f "$cfg/annotations.tsv" ] || continue
  echo ">>> $(basename "$(dirname "$cfg")")/$(basename "$cfg")  $(date '+%H:%M:%S')"
  conda run --no-capture-output -n "${INFERCNV_ENV:-infercnv-r}" Rscript \
    "$SCRIPT_DIR/sensitivity_reference.R" --config-dir "${cfg%/}" --threads "$THREADS" \
    || echo "[FAIL] ${cfg%/}" >&2
done

conda run --no-capture-output -n "${PREP_ENV:-benchmark-py-r}" \
  python3 "$SCRIPT_DIR/sensitivity_reference.py" --mode analyze
