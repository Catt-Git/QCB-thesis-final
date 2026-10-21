#!/usr/bin/env bash
#SBATCH --job-name=integration
#SBATCH --output=logs/integration_%A_%a.out
#SBATCH --error=logs/integration_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=1-17%3
#
# 02_2 integration: SLURM array (CPU, GPU was not present in my case), one task per grid row, at most 3 at once.
# DRVI is not in the array (see submit_drvi.slurm).
#
#   INTEGRATION_ENV  every method except scGen   [default: benchmark-py-r]
#   SCGEN_ENV        scGen only                  [default: scgen-py; or run scGen locally]
#
# Usage:
#   cd 02_integration_benchmark/02_2_integration   # submit FROM this folder
#   export DATA_DIR=<data root>
#   bash run_all.sh --slurm                        # usual route: a filtered subset, same %3 throttle
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_integration.slurm   # direct use, the whole array

set -euo pipefail
mkdir -p logs

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
INTEGRATION_ENV="${INTEGRATION_ENV:-benchmark-py-r}"
SCGEN_ENV="${SCGEN_ENV:-scgen-py}"

# SLURM runs a spool copy, so locate the repo via INTEGRATION_DIR or the submit dir.
SCRIPT_DIR="${INTEGRATION_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
[ -f "$GRID" ] || { echo "grid not found: $GRID (set INTEGRATION_DIR, or sbatch from 02_2_integration/)" >&2; exit 1; }
HVG_RDS="$DATA_DIR/shiao_hvg_2k_unintegrated_list.rds"

has_embed() { case ",$1," in *,embed,*) return 0 ;; *) return 1 ;; esac }

# Checkpoint dir for trained methods, as in run_all.sh; empty otherwise.
set_model_arg() {
  model_arg=()
  case "$1" in
    scvi|scanvi|scgen) model_arg=(--model-dir "$DATA_DIR/02_${1}") ;;
  esac
}

# Pick this task's grid row (index = row number, matching run_all.sh's expansion).
JOB=$(awk -F'\t' 'NR>1 && $1!="" { idx++; if (idx == '"$SLURM_ARRAY_TASK_ID"') print }' "$GRID")
[ -n "$JOB" ] || { echo "no grid row for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r run_id method language env scaling input output types ref_h5ad hvgs <<< "$JOB"

# Reference batches for Seurat CCA/RPCA (not the grid's `reference` column).
SEURAT_REF_ARG=()
if [ -n "${SEURAT_REFERENCE:-}" ]; then SEURAT_REF_ARG=(--reference "$SEURAT_REFERENCE"); fi

if [ "$language" = "notebook" ]; then
  echo "[skip] $run_id is a notebook (DRVI); submit it with submit_drvi.slurm"; exit 0
fi

# scGen has its own env; everything else uses INTEGRATION_ENV.
if [ "$method" = "scgen" ]; then ENVUSE="$SCGEN_ENV"; else ENVUSE="$INTEGRATION_ENV"; fi

in_path="$DATA_DIR/$input"
out_path="$DATA_DIR/$output"
emb_path="$DATA_DIR/02_embeddings/${run_id}.npy"
mkdir -p "$(dirname "$out_path")" "$(dirname "$emb_path")"

echo ">>> [$SLURM_ARRAY_TASK_ID] $run_id  ($method, $language, env=$ENVUSE)"

# set -u off for the conda block: base's openjdk activate script reads unset JAVA_HOME.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$ENVUSE"
set -u

case "$language" in
  python)
    set_model_arg "$method"
    if [ "$method" = "scgen" ]; then
      python "$SCRIPT_DIR/run_scgen.py" -i "$in_path" -o "$out_path" "${model_arg[@]}"
    else
      emb_arg=(); has_embed "$types" && emb_arg=(--emb-out "$emb_path")
      python "$SCRIPT_DIR/run_integration.py" -m "$method" -i "$in_path" -o "$out_path" \
        "${emb_arg[@]}" "${model_arg[@]}"
    fi
    ;;
  R)
    rds_out="$DATA_DIR/02_integration/${run_id}.rds"
    Rscript "$SCRIPT_DIR/run_integration.R" -m "$method" -i "$in_path" -o "$rds_out" \
      -v "$HVG_RDS" "${SEURAT_REF_ARG[@]}"
    emb_arg=(); has_embed "$types" && emb_arg=(--emb-out "$emb_path")
    python "$SCRIPT_DIR/rds_to_h5ad.py" -i "$rds_out" -o "$out_path" --types "$types" "${emb_arg[@]}"
    # Drop the large .rds once converted (KEEP_RDS=1 to keep it).
    if [ "${KEEP_RDS:-0}" = "1" ]; then echo "[keep] $rds_out"; else
      echo "[clean] removing intermediate $rds_out"; rm -f "$rds_out"; fi
    ;;
  *)
    echo "[error] unknown language '$language' for $run_id" >&2; exit 1 ;;
esac

echo "[ok] $run_id -> $out_path"
