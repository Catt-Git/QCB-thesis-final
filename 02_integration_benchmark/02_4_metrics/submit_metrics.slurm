#!/usr/bin/env bash
#SBATCH --job-name=scib_metrics
#SBATCH --output=logs/metrics_%A_%a.out
#SBATCH --error=logs/metrics_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --array=1-22%5
#
# 02_4 metrics: SLURM array, one task per (run, output type) (18 rows -> 22 tasks),
# 5 at a time. 96G is enough except for the bbknn knn tasks (leiden on 8.7e8
# nonzeros peaked at 197 GiB): submit those with more memory (last usage line).
#
# Usage:
#   cd 02_integration_benchmark/02_4_metrics               # submit FROM this folder
#   export DATA_DIR=<data root>
#   bash run_all_metrics.sh --slurm                        # usual route
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_metrics.slurm   # direct use, the whole array
#   sbatch --mem=300G --array=1 ... submit_metrics.slurm   # bbknn_unscaled knn

set -euo pipefail
mkdir -p logs

ENV="${METRICS_ENV:-benchmark-py-r}"
TASK="shiao"
BATCH_KEY="cohort"; LABEL_KEY="cell_type"; ORGANISM="human"
# Leiden flavor for NMI/ARI (see metrics_shared.use_igraph_leiden).
CLUSTER_FLAVOR="${CLUSTER_FLAVOR:-igraph}"

# SLURM runs a spool copy, so locate the repo via METRICS_DIR or the submit dir.
SCRIPT_DIR="${METRICS_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
[ -f "$GRID" ] || { echo "grid not found: $GRID (set METRICS_DIR, or sbatch from 02_4_metrics/)" >&2; exit 1; }

# DATA_DIR from --export (cluster_env.sh).
. "$SCRIPT_DIR/../utils/cluster_env.sh"
resolve_data_dir || exit 1
ROOT="$DATA_DIR/02_metrics"

# Grid expanded to one line per (run, type); line SLURM_ARRAY_TASK_ID is this task.
JOB=$(awk -F'\t' 'NR>1 && $1!="" {
        n=split($8, ts, ",");
        for (i=1;i<=n;i++) print $1"\t"$2"\t"$5"\t"$7"\t"ts[i]"\t"$9"\t"$10
      }' "$GRID" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$JOB" ] || { echo "no job for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r run_id method scaling output type reference hvgs <<< "$JOB"

int_path="$DATA_DIR/$output"
ref_path="$DATA_DIR/$reference"
out_dir="$ROOT/$TASK/metrics/$scaling/hvg"
csv="$out_dir/${method}_${type}.csv"
mkdir -p "$out_dir"

echo ">>> [$SLURM_ARRAY_TASK_ID] $run_id  ($method, $type)"

# Unbuffered, so an OOM kill does not lose the progress log.
export PYTHONUNBUFFERED=1

# LISI writes tens of GB to TMPDIR on BBKNN; it must be disk-backed (METRICS_TMPDIR).
export TMPDIR="${METRICS_TMPDIR:-/tmp}/scib_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "[node] $(hostname)  mem=${SLURM_MEM_PER_NODE:-?}M  cpus=${SLURM_CPUS_PER_TASK:-?}  TMPDIR=$TMPDIR"

# set -u off for the conda block: base's openjdk activate script reads unset JAVA_HOME.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$ENV"
set -u

python "$SCRIPT_DIR/check_integrations.py" -i "$int_path" -u "$ref_path" \
  --types "$type" -b "$BATCH_KEY" -l "$LABEL_KEY"

python "$SCRIPT_DIR/metrics.py" -u "$ref_path" -i "$int_path" -o "$csv" \
  -m "$method" --type "$type" --scaling "$scaling" -b "$BATCH_KEY" -l "$LABEL_KEY" \
  --organism "$ORGANISM" --hvgs "$hvgs" --cluster-flavor "$CLUSTER_FLAVOR" -v

echo "[ok] $run_id $type -> $csv"
