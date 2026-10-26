#!/usr/bin/env bash
#
# 03 utils: --slurm for the phase drivers. Sourced, never executed:
#   . "$PHASE_DIR/utils/slurm.sh"
#
# submit_driver turns a whole driver into one job (utils/submit_driver.slurm) that re-runs
# it without --slurm, so resuming works on the cluster as it does locally.
#
#   JOB_PARTITION   partition                          [default: long]
#   JOB_CPUS        overrides the driver's CPU count
#   JOB_MEM         overrides the driver's memory
#   JOB_ENV         conda env activated in the job     [default: benchmark-py-r]

SLURM_UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# submit_driver <job name> <cpus> <mem> <log dir> <dry run 0|1> <driver> [driver args...]
submit_driver() {
  local name="$1" cpus="${JOB_CPUS:-$2}" mem="${JOB_MEM:-$3}" log_dir="$4" dry="$5" driver="$6"
  shift 6
  local args=() a
  for a in "$@"; do [ "$a" = "--slurm" ] || args+=("$a"); done

  local partition="${JOB_PARTITION:-long}"
  # Absolute log paths, so logs land in logs/ wherever sbatch is called from.
  local opts=(--job-name="$name" --partition="$partition" --nodes=1
              --cpus-per-task="$cpus" --mem="$mem" --chdir="$(dirname "$driver")"
              --output="$log_dir/${name}_%j.out" --error="$log_dir/${name}_%j.err"
              --export="ALL,DATA_DIR=$DATA_DIR")
  local wrapper="$SLURM_UTILS_DIR/submit_driver.slurm"

  if [ "$dry" -eq 1 ]; then
    echo "[dry] sbatch ${opts[*]} $wrapper $driver ${args[*]:-}"
    return 0
  fi
  # logs/ must exist before submission: SLURM opens --output before the job starts.
  mkdir -p "$log_dir"
  local j
  # "${a[@]+...}": an empty array passes no argument.
  j="$(sbatch --parsable "${opts[@]}" "$wrapper" "$driver" ${args[@]+"${args[@]}"})"
  echo "submitted $(basename "$driver") as job $j ($partition, $cpus CPUs, $mem, env ${JOB_ENV:-benchmark-py-r})"
  echo "logs : $log_dir/${name}_$j.out (and .err)"
  echo "watch: squeue -j $j ; tail -f $log_dir/${name}_$j.out"
}
