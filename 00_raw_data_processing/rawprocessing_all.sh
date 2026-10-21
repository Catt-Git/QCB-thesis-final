#!/usr/bin/env bash
#
# 00 raw data processing: the whole phase in one command, on the cluster.
#
# Submits the seven SLURM jobs of the phase as a chain, each one starting only when
# the previous one has succeeded (--dependency=afterok):
#
#   1. metadata    00_1_metadata/submit_metadata.slurm          -> ena_full_metadata.tsv (repo)
#   2. samplemap   00_2_sample_mapping/submit_sample_map.slurm  -> sample_map_gex_final.tsv (repo)
#   3. download    00_3_download/submit_download.slurm          -> $DATA_DIR/fastq_raw/ (array)
#   4. verify      00_3_download/submit_verify_integrity.slurm  -> verify_integrity_<jobid>.{out,report}
#   5. rename      00_4_rename/submit_rename.slurm              -> $DATA_DIR/renamed_fastq/
#   6. cellranger  00_5_cellranger/submit_cellranger.slurm      -> $DATA_DIR/cellranger_out/
#   7. h5ad        00_6_build_h5ad/submit_h5ad_concat.slurm     -> $DATA_DIR/all_samples_combined.h5ad
#
# Resuming is the default: the [have]/[run] decision is taken here, at submit time,
# from what is on disk:
#
#   metadata, samplemap  the TSV exists
#   download             only the missing runs are submitted, as --array=<their lines>
#   verify               the latest verify_integrity_<jobid>.out says "[ok] no corrupted
#                        file" and no FASTQ is newer than it
#   rename               every run has both renamed copies, none older than its source
#   cellranger           every sample in unique_samples.txt has its filtered matrix
#   h5ad                 all_samples_combined.h5ad exists
#
# Once all_samples_combined.h5ad exists the phase is complete and nothing before it is
# looked at, so deleting the FASTQ to free space (README) does not trigger a 2.4 TB
# re-download. Delete the .h5ad, name the steps, or pass --force to run steps again.
#
# A step that writes data (download, rename, cellranger) drags every later step along
# once it is submitted, whatever is on disk: re-downloaded FASTQ have to be checked and
# renamed again, and new matrices mean a new .h5ad. --force submits every step, but
# the download still fetches only what is missing: re-fetching a run means deleting
# its FASTQ first (README, "Repairing a single run").
#
# Usage:
#   export DATA_DIR=<data root>
#   bash rawprocessing_all.sh                    # submit what is missing, chained
#   bash rawprocessing_all.sh --dry-run          # print the plan and the sbatch commands, submit nothing
#   bash rawprocessing_all.sh --force            # submit every step
#   bash rawprocessing_all.sh download verify    # only the named step(s), chained in phase order
#
# Step names: metadata, samplemap, download, verify, rename, cellranger, h5ad.
# REF and DOWNLOAD_ENV reach the jobs when exported (sbatch --export=ALL); otherwise
# the .slurm defaults apply. Each job logs next to its own .slurm, as when submitted by
# hand; the plan and the job ids also go to 00_raw_data_processing/logs/rawprocessing_all_<timestamp>.log.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
# Concurrent download tasks: the same %4 as the #SBATCH --array line, which the
# --array passed here replaces.
THROTTLE=4

# One row per step: name, folder it is submitted from (the .slurm files anchor on
# SLURM_SUBMIT_DIR), batch script.
STEP_NAMES=(metadata samplemap download verify rename cellranger h5ad)
STEP_DIRS=(
  00_1_metadata
  00_2_sample_mapping
  00_3_download
  00_3_download
  00_4_rename
  00_5_cellranger
  00_6_build_h5ad
)
STEP_SLURMS=(
  submit_metadata.slurm
  submit_sample_map.slurm
  submit_download.slurm
  submit_verify_integrity.slurm
  submit_rename.slurm
  submit_cellranger.slurm
  submit_h5ad_concat.slurm
)
DATA_WRITERS=(download rename cellranger)

META="$SCRIPT_DIR/00_1_metadata/ena_full_metadata.tsv"
MAP="$SCRIPT_DIR/00_2_sample_mapping/sample_map_gex_final.tsv"

FORCE=0; DRY_RUN=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,47p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the data root on the cluster (needs ~4.8TB free for fastq_raw/ + renamed_fastq/)}"
FASTQ_DIR="$DATA_DIR/fastq_raw"
RENAMED_DIR="$DATA_DIR/renamed_fastq"
CR_DIR="$DATA_DIR/cellranger_out"
H5AD="$DATA_DIR/all_samples_combined.h5ad"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Validate the step names before looking at the runs on disk, so a typo fails at once.
for requested in "${STEP_FILTER[@]:-}"; do
  [ -n "$requested" ] || continue
  in_list "$requested" "${STEP_NAMES[@]}" || {
    echo "unknown step '$requested'; valid steps: ${STEP_NAMES[*]}" >&2; exit 1; }
done

scheduled() {  # is this step name wanted on this invocation?
  local name="$1"
  [ ${#STEP_FILTER[@]} -eq 0 ] && return 0
  in_list "$name" "${STEP_FILTER[@]}"
}

# srr <TAB> sample_id, one row per line of the map. Read through awk rather than with
# `read` on the raw map: column 4 (sub-number) is usually empty, and `read` collapses
# consecutive tabs, which would shift sample_id out of column 6.
map_rows() { awk -F'\t' '{print $1"\t"$6}' "$MAP"; }

# Line numbers of the map whose run is not fully downloaded: line N is array task N.
missing_lines() {
  local n=0 srr sample
  while IFS=$'\t' read -r srr sample; do
    n=$((n + 1))
    [ -n "$srr" ] || continue
    [ -f "$FASTQ_DIR/${srr}_1.fastq.gz" ] && [ -f "$FASTQ_DIR/${srr}_2.fastq.gz" ] && continue
    echo "$n"
  done < <(map_rows)
}

to_ranges() {  # one number per line -> 1-3,7,9-10
  awk 'NR == 1 { s = p = $1; next }
       $1 == p + 1 { p = $1; next }
       { out = out (out ? "," : "") (s == p ? s : s "-" p); s = p = $1 }
       END { if (NR) print out (out ? "," : "") (s == p ? s : s "-" p) }'
}

download_complete() { [ -s "$MAP" ] && [ -z "$(missing_lines)" ]; }

# The latest verify job passed, and nothing was downloaded after it.
verify_done() {
  local out
  out="$(ls -t "$SCRIPT_DIR"/00_3_download/verify_integrity_*.out 2>/dev/null | head -1 || true)"
  [ -n "$out" ] && grep -q '^\[ok\] no corrupted file' "$out" || return 1
  [ -z "$(find "$FASTQ_DIR" -maxdepth 1 -name '*.fastq.gz' -newer "$out" -print -quit 2>/dev/null)" ]
}

# rsync -a keeps the source mtime, so a source newer than its copy was re-downloaded
# after the renaming and the copy is stale.
rename_done() {
  [ -s "$RENAMED_DIR/unique_samples.txt" ] && [ -s "$MAP" ] || return 1
  local srr sample r src dst
  while IFS=$'\t' read -r srr sample; do
    [ -n "$srr" ] || continue
    for r in 1 2; do
      src="$FASTQ_DIR/${srr}_${r}.fastq.gz"
      dst="$RENAMED_DIR/${sample}_S1_L001_R${r}_001.fastq.gz"
      [ -f "$dst" ] || return 1
      if [ "$src" -nt "$dst" ]; then return 1; fi
    done
  done < <(map_rows)
}

cellranger_done() {
  [ -s "$RENAMED_DIR/unique_samples.txt" ] || return 1
  bash "$SCRIPT_DIR/utils/check_cellranger_completeness.sh" \
    "$RENAMED_DIR/unique_samples.txt" "$CR_DIR" >/dev/null 2>&1
}

step_done() {  # step index -> 0 if its result is already on disk
  [ "$FORCE" -eq 1 ] && return 1
  case "${STEP_NAMES[$1]}" in
    metadata)   [ -s "$META" ] ;;
    samplemap)  [ -s "$MAP" ] ;;
    download)   download_complete ;;
    verify)     verify_done ;;
    rename)     rename_done ;;
    cellranger) cellranger_done ;;
    h5ad)       [ -f "$H5AD" ] ;;
  esac
}

input_ok() {  # step index -> 0 if what it reads is already on disk
  case "${STEP_NAMES[$1]}" in
    metadata)      return 0 ;;
    samplemap)     [ -s "$META" ] ;;
    download)      [ -s "$MAP" ] ;;
    verify|rename) download_complete ;;
    cellranger)    [ -s "$RENAMED_DIR/unique_samples.txt" ] ;;
    h5ad)          cellranger_done ;;
  esac
}

input_hint() {
  case "${STEP_NAMES[$1]}" in
    samplemap)     echo "missing $META; add the metadata step" ;;
    download)      echo "missing $MAP; add the samplemap step" ;;
    verify|rename) echo "runs missing from $FASTQ_DIR; add the download step (utils/check_download_completeness.sh lists them)" ;;
    cellranger)    echo "missing $RENAMED_DIR/unique_samples.txt; add the rename step" ;;
    h5ad)          echo "Cell Ranger matrices missing under $CR_DIR; add the cellranger step (utils/check_cellranger_completeness.sh lists them)" ;;
  esac
}

# Plan: one of run / have / skip per step. A step's input may be missing only if the
# step right before it in the phase is about to produce it; otherwise fail now rather
# than after a day in the queue. In --dry-run it only warns, so the plan stays
# inspectable.
PLAN=()
DL_ARRAY=""        # --array value for the download, empty = the #SBATCH default
drag=0; prev_run=0; missing_input=0; warnings=()
PHASE_DONE=0
[ -f "$H5AD" ] && [ "$FORCE" -eq 0 ] && [ ${#STEP_FILTER[@]} -eq 0 ] && PHASE_DONE=1
for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  if ! scheduled "$name"; then PLAN[i]=skip; prev_run=0; continue; fi
  if [ "$PHASE_DONE" -eq 1 ]; then PLAN[i]=have; continue; fi

  if [ "$drag" -eq 1 ] || ! step_done "$i"; then PLAN[i]=run; else PLAN[i]=have; fi

  if [ "${PLAN[$i]}" = run ] && [ "$name" = download ]; then
    if [ -s "$MAP" ]; then
      lines="$(missing_lines)"
      if [ -z "$lines" ]; then
        PLAN[i]=have     # forced or not, there is nothing left to fetch
      else
        DL_ARRAY="$(to_ranges <<< "$lines")%$THROTTLE"
        # A run fetched again after Cell Ranger already saw it: step 5 skips a sample
        # whose filtered matrix exists, so the old matrix would survive.
        while IFS=$'\t' read -r srr sample; do
          [ -d "$CR_DIR/$sample" ] && \
            warnings+=("$srr is re-downloaded but $CR_DIR/$sample exists: delete it, or step 5 keeps the old matrix")
        done < <(map_rows | awk -v want="$(tr '\n' ' ' <<< "$lines")" \
                   'BEGIN { n = split(want, a, " "); for (k = 1; k <= n; k++) w[a[k]] = 1 } NR in w')
      fi
    else
      warnings+=("no sample map yet: the download is submitted with the #SBATCH --array of submit_download.slurm")
    fi
  fi

  if [ "${PLAN[$i]}" = run ]; then
    in_list "$name" "${DATA_WRITERS[@]}" && drag=1
    if [ "$prev_run" -eq 0 ] && ! input_ok "$i"; then
      echo "[error] $name: $(input_hint "$i")" >&2
      missing_input=1
    fi
    prev_run=1
  else
    prev_run=0
  fi
done
if [ "$missing_input" -eq 1 ]; then
  [ "$DRY_RUN" -eq 0 ] && exit 1
  echo "[warn] --dry-run: continuing anyway, the plan below would fail as it stands"
fi

if [ "$DRY_RUN" -eq 0 ]; then
  command -v sbatch >/dev/null || { echo "[error] sbatch not found: run this on the cluster's login node" >&2; exit 1; }
  mkdir -p "$LOG_DIR"
  LOG_FILE="$LOG_DIR/rawprocessing_all_$(date +%Y%m%d_%H%M%S).log"
  exec > >(tee -a "$LOG_FILE") 2>&1
fi

echo "DATA_DIR : $DATA_DIR"
[ "$DRY_RUN" -eq 0 ] && echo "log      : $LOG_FILE"
echo
[ "$PHASE_DONE" -eq 1 ] && \
  echo "[have] $H5AD exists: the phase is complete (delete it, or name the steps, to run them again)"
for w in "${warnings[@]+"${warnings[@]}"}"; do echo "[warn] $w"; done

# Submit in phase order, each job waiting on the last one submitted. kill-on-invalid-dep:
# when a job fails, the ones waiting on it are cancelled instead of sitting in the queue
# as DependencyNeverSatisfied.
JOBS=(); prev_job=""; n_run=0; n_have=0; n_skip=0
for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  case "${PLAN[$i]}" in
    skip) n_skip=$((n_skip + 1)); continue ;;
    have)
      if [ "$PHASE_DONE" -eq 1 ]; then
        :     # one line above says it all; the per-step checks were not even run
      elif [ "$name" = download ]; then
        echo "[have] download: every run is on disk (to fetch one again, delete its FASTQ first, see README)"
      else
        echo "[have] $name: already done, not submitted"
      fi
      n_have=$((n_have + 1)); continue ;;
  esac

  args=(--parsable --export="ALL,DATA_DIR=$DATA_DIR")
  [ -n "$prev_job" ] && args+=(--dependency="afterok:$prev_job" --kill-on-invalid-dep=yes)
  [ "$name" = download ] && [ -n "$DL_ARRAY" ] && args+=(--array="$DL_ARRAY")

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] $name: (cd ${STEP_DIRS[$i]} && sbatch ${args[*]} ${STEP_SLURMS[$i]})"
    prev_job="<$name>"
  else
    j="$(cd "$SCRIPT_DIR/${STEP_DIRS[$i]}" && sbatch "${args[@]}" "${STEP_SLURMS[$i]}")" || {
      echo "[error] sbatch failed for $name; already queued: ${JOBS[*]:-none} (scancel them or let them run)" >&2
      exit 1; }
    j="${j%%;*}"     # --parsable prints "jobid[;cluster]"
    note=""
    [ "$name" = download ] && [ -n "$DL_ARRAY" ] && note=", array $DL_ARRAY"
    echo "[sub] $name: job $j${prev_job:+ (after $prev_job)}$note"
    prev_job="$j"; JOBS+=("$j")
  fi
  n_run=$((n_run + 1))
done

echo "==================================================================="
echo "$n_run step(s) submitted, $n_have already done, $n_skip filtered out"
[ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
  echo "      (--force to submit the $n_have done step(s) again)"
if [ ${#JOBS[@]} -gt 0 ]; then
  ids="$(IFS=,; echo "${JOBS[*]}")"
  echo "watch: squeue -j $ids ; sacct -j $ids --format=JobID,JobName,State,Elapsed,ExitCode"
fi
[ "${PLAN[6]}" = run ] && echo "next : 01_pre_processing on $H5AD"
exit 0
