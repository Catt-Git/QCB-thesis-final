#!/usr/bin/env bash
#
# 03_1 sensitivity_resolution: does the per-patient verdict depend on the subcluster resolution?
#
# Runs run_infercnv.R at leiden 0.005 / 0.01 / 0.02 / 0.05 on four cohorts: Patient16
# (tumour-rich), Patient43 (intermediate), Patient64 (small stromal null), Patient52 (no
# detectable tumour). The existing summaries are restored after each run.
#
# Output: $DATA_DIR/03_tum/sensitivity/resolution/<cohort>_res<r>.csv (+ .log)
#
# Usage:
#   export DATA_DIR=<data root>
#   bash sensitivity_resolution.sh                # 16 runs, ~75 min, resumes
#   bash sensitivity_resolution.sh --threads 8
#   python3 sensitivity_resolution.py             # then the tables

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${DATA_DIR:?set DATA_DIR}"

THREADS="${THREADS:-20}"
COHORTS=(Patient64 Patient43 Patient16 Patient52)
RESOLUTIONS=(0.005 0.01 0.02 0.05)
[ "${1:-}" = "--threads" ] && THREADS="$2"

OUT="$DATA_DIR/03_tum/sensitivity/resolution"
SUMMARY="$DATA_DIR/03_tum/summary"
mkdir -p "$OUT"

for res in "${RESOLUTIONS[@]}"; do
  for coh in "${COHORTS[@]}"; do
    dest="$OUT/${coh}_res${res}.csv"
    [ -f "$dest" ] && { echo "[have] $coh res=$res"; continue; }
    # run_infercnv.R writes to summary/: move the existing summary aside and restore it after.
    keep=""
    if [ -f "$SUMMARY/${coh}_cnv.csv" ]; then keep="$SUMMARY/${coh}_cnv.csv.sensitivity_bak"
      mv "$SUMMARY/${coh}_cnv.csv" "$keep"; fi
    echo ">>> $coh res=$res  $(date '+%H:%M:%S')"
    if conda run --no-capture-output -n "${INFERCNV_ENV:-infercnv-r}" Rscript \
         "$SCRIPT_DIR/run_infercnv.R" --cohort "$coh" --threads "$THREADS" --force \
         --leiden-resolution "$res" > "$OUT/${coh}_res${res}.log" 2>&1 \
       && [ -f "$SUMMARY/${coh}_cnv.csv" ]; then
      mv "$SUMMARY/${coh}_cnv.csv" "$dest"; echo "[ok] $coh res=$res"
    else
      echo "[FAIL] $coh res=$res - see $OUT/${coh}_res${res}.log" >&2
      tail -5 "$OUT/${coh}_res${res}.log" >&2
    fi
    [ -n "$keep" ] && mv "$keep" "$SUMMARY/${coh}_cnv.csv"
  done
done
echo "done - now: python3 sensitivity_resolution.py"
