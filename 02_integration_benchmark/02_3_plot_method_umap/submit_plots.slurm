#!/usr/bin/env bash
#SBATCH --job-name=scib_plots
#SBATCH --output=logs/plots_%A_%a.out
#SBATCH --error=logs/plots_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --array=1-17%5
#
# 02_3 method plots: SLURM array, one task per grid row (index = row number),
# 5 at a time. 96G is enough except for the two bbknn rows (see last usage line).
#
# Usage:
#   cd 02_integration_benchmark/02_3_plot_method_umap      # submit FROM this folder
#   export DATA_DIR=<data root>
#   bash plot_all.sh --slurm                               # usual route
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_plots.slurm   # direct use, the whole array
#   sbatch --mem=300G --array=1,2 ... submit_plots.slurm   # the bbknn rows

set -euo pipefail
mkdir -p logs

ENV="${PLOT_ENV:-benchmark-py-r}"

# SLURM runs a spool copy, so locate the repo via PLOT_DIR or the submit dir.
SCRIPT_DIR="${PLOT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
PHASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GRID="$PHASE_DIR/benchmark_grid.tsv"
FIG_DIR="$PHASE_DIR/figures"
PLOTTER="$SCRIPT_DIR/plot_methods_umaps.py"
[ -f "$GRID" ] || { echo "grid not found: $GRID (set PLOT_DIR, or sbatch from 02_3_plot_method_umap/)" >&2; exit 1; }

# DATA_DIR from --export (cluster_env.sh).
. "$SCRIPT_DIR/../utils/cluster_env.sh"
resolve_data_dir || exit 1

# The SLURM_ARRAY_TASK_ID-th data row of the grid is this task.
# Columns: run_id method scaling output types reference.
JOB=$(awk -F'\t' 'NR>1 && $1!="" { print $1"\t"$2"\t"$5"\t"$7"\t"$8"\t"$9 }' "$GRID" \
      | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$JOB" ] || { echo "no job for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r run_id method scaling output types reference <<< "$JOB"

# Must match plot_type() in plot_all.sh.
case ",$types," in
  *,embed,*) ptype=embed ;;
  *,full,*)  ptype=full ;;
  *,knn,*)   ptype=knn ;;
  *)         ptype="$types" ;;
esac

int_path="$DATA_DIR/$output"
ref_path="$DATA_DIR/$reference"
out_dir="$FIG_DIR/$run_id"

echo ">>> [$SLURM_ARRAY_TASK_ID] $run_id  ($method, --type $ptype)"
echo "[node] $(hostname)  mem=${SLURM_MEM_PER_NODE:-?}M  cpus=${SLURM_CPUS_PER_TASK:-?}"
echo "[out ] $out_dir"

# Unbuffered, so an OOM kill does not lose the progress log.
export PYTHONUNBUFFERED=1
# Headless node.
export MPLBACKEND=Agg

# set -u off for the conda block: base's openjdk activate script reads unset JAVA_HOME.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$ENV"
set -u

# --recompute-umap forwarded by plot_all.sh via --export.
UMAP_ARGS=()
[ "${PLOT_RECOMPUTE_UMAP:-0}" = "1" ] && UMAP_ARGS+=(--recompute-umap)

python "$PLOTTER" -m "$method" --type "$ptype" \
  -i "$int_path" -u "$ref_path" -o "$out_dir" ${UMAP_ARGS[@]+"${UMAP_ARGS[@]}"}

echo "[ok] $run_id -> $out_dir"
