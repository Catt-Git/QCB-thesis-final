#!/usr/bin/env bash
#
# 02 utils: local smoke test of the 02_4 pipeline (check_integrations.py,
# metrics.py, merge_metrics.py) on the small 02_2 outputs in $SMOKE_DIR/smoke_out/,
# against smoke_hvg.h5ad. Every expected metric must be finite. Unscaled rows only;
# DRVI is skipped. SMOKE_DIR defaults to utils/smoke_test/.
#
# Usage:
#   conda activate benchmark-py-r                          # so rpy2 can find R
#   bash smoke_test_metrics_pipeline.sh --dry-run          # print what would run, do nothing
#   bash smoke_test_metrics_pipeline.sh
#   bash smoke_test_metrics_pipeline.sh --method harmony   # one method only
#   bash smoke_test_metrics_pipeline.sh bbknn_unscaled     # one run_id only
#
# Outputs (all disposable): $SMOKE_DIR/smoke_metrics/ + smoke_metrics_merged.csv.

set -euo pipefail

UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SMOKE_DIR="${SMOKE_DIR:-$UTILS_DIR/smoke_test}"
METRICS_DIR="$UTILS_DIR/../02_4_metrics"
GRID="$UTILS_DIR/../benchmark_grid.tsv"
ENV="${METRICS_ENV:-benchmark-py-r}"
TASK="shiao"
BATCH_KEY="cohort"; LABEL_KEY="cell_type"; ORGANISM="human"

ROOT="$SMOKE_DIR/smoke_metrics"
MERGED="$SMOKE_DIR/smoke_metrics_merged.csv"

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

conda_guarded() { set +u; conda "$@"; set -u; }

# Parse the command line: --method, run_id filters, --no-check, --dry-run.
NO_CHECK=0; DRY_RUN=0
METHOD_FILTER=(); FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --no-check) NO_CHECK=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --method) METHOD_FILTER+=("${2:-}"); shift 2 ;;
    --method=*) METHOD_FILTER+=("${1#*=}"); shift ;;
    *) FILTER+=("$1"); shift ;;
  esac
done

# Rows matching every filter (method, run_id).
want() {
  local run_id="$1" method="$2" f matched
  if [ ${#METHOD_FILTER[@]} -gt 0 ]; then
    matched=0
    for f in "${METHOD_FILTER[@]}"; do [ "$f" = "$method" ] && matched=1; done
    [ "$matched" -eq 1 ] || return 1
  fi
  [ ${#FILTER[@]} -eq 0 ] && return 0
  for f in "${FILTER[@]}"; do [ "$f" = "$run_id" ] && return 0; done
  return 1
}

if [ "$DRY_RUN" -eq 0 ]; then
  set +u; eval "$(conda shell.bash hook)"; set -u
  conda_guarded activate "$ENV"
fi

n_done=0; n_skip=0

{
  read -r _header
  while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
    [ -n "$run_id" ] || continue
    [ "$scaling" = "unscaled" ] || { n_skip=$((n_skip + 1)); continue; }
    if ! want "$run_id" "$method"; then n_skip=$((n_skip + 1)); continue; fi
    if [ "$language" = "notebook" ]; then
      echo "[skip] $run_id: notebook method, no smoke output"; n_skip=$((n_skip + 1)); continue
    fi

    int_path="$SMOKE_DIR/smoke_out/${method}.h5ad"
    ref_path="$SMOKE_DIR/smoke_hvg.h5ad"
    if [ ! -f "$int_path" ]; then
      echo "[skip] $run_id: no smoke output ($int_path)"; n_skip=$((n_skip + 1)); continue
    fi

    out_dir="$ROOT/$TASK/metrics/$scaling/hvg"

    if [ "$DRY_RUN" -eq 1 ]; then
      echo "[dry] $run_id  types=$types"
      IFS=',' read -ra tps <<< "$types"
      for t in "${tps[@]}"; do echo "      $method $t -> $out_dir/${method}_${t}.csv"; done
      n_done=$((n_done + 1)); continue
    fi

    mkdir -p "$out_dir"
    echo "==================================================================="
    echo ">>> $run_id  ($method, types=$types)"
    echo "==================================================================="

    if [ "$NO_CHECK" -eq 0 ]; then
      python "$METRICS_DIR/check_integrations.py" \
        -i "$int_path" -u "$ref_path" --types "$types" -b "$BATCH_KEY" -l "$LABEL_KEY"
    fi

    IFS=',' read -ra tps <<< "$types"
    for t in "${tps[@]}"; do
      t="$(echo "$t" | tr -d '[:space:]')"
      csv="$out_dir/${method}_${t}.csv"
      echo "--- metrics: $method $t ---"
      python "$METRICS_DIR/metrics.py" \
        -u "$ref_path" -i "$int_path" -o "$csv" \
        -m "$method" --type "$t" --scaling "$scaling" -b "$BATCH_KEY" -l "$LABEL_KEY" \
        --organism "$ORGANISM" --hvgs "$hvgs" -v
    done

    n_done=$((n_done + 1))
    echo "[ok] $run_id"
  done
} < "$GRID"

if [ "$DRY_RUN" -eq 0 ]; then
  echo "==================================================================="
  echo ">>> merging per-run CSVs"
  python "$METRICS_DIR/merge_metrics.py" \
    -o "$MERGED" -r "$ROOT" --glob "$ROOT/$TASK/metrics/*/hvg/*.csv" || \
    echo "[warn] merge produced nothing (no CSVs yet?)"

  # Fail if a metric expected for the run's output type has no finite value.
  if [ -f "$MERGED" ]; then
    echo ">>> checking that every expected metric is finite"
    python - "$MERGED" "$UTILS_DIR" <<'PY'
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, sys.argv[2])
from metrics_shared import EXPECTED

table = pd.read_csv(sys.argv[1], index_col=0)
bad = []
for run, row in table.iterrows():
    expected = list(EXPECTED[run.rsplit("_", 1)[1]])
    if "/scaled/" in run and "PCR_batch" in expected:
        expected.remove("PCR_batch")
    bad += [f"{run}: {m}" for m in expected if not np.isfinite(row.get(m, np.nan))]
if bad:
    print("[fail] metrics without a value:\n  " + "\n  ".join(bad))
    sys.exit(1)
print(f"[ok] {len(table)} run(s), every expected metric is finite")
PY
  fi
  conda_guarded deactivate
fi

echo "==================================================================="
echo "done: $n_done run(s) scored, $n_skip skipped   ->  $MERGED"
