#!/usr/bin/env python3

# build_requested_metadata.py
# Build sample_metadata_final.tsv from sample_map_gex_final.tsv.
# Output columns (join key downstream = sample_id): 
# - srr (SRA/ENA run accession)
# - sample_id (Cell Ranger folder, e.g. P01_A_P / P03_B2_P)
# - batch (original library label / batch_key, e.g. h01A_P)
# - cohort (PatientNN, zero-padded)
# - treatment (mapped from timepoint via TREATMENT_MAP)
# - fraction (CD45+ for P / CD45- for N)
# - dataset_origin (immune for P / non_immune for N)
# - response (R1/R2/NR, PATIENT-LEVEL, from RESPONSE_MAP).

# Usage:
#   python3 build_requested_metadata.py \
#       --sample-map ../00_2_sample_mapping/sample_map_gex_final.tsv \
#       --out        sample_metadata_final.tsv

# On the cluster (how the thesis run was done) this is not launched on its own: submit_h5ad_concat.slurm
# runs it first and then build_combined_h5ad.py, so the merge always reads a metadata table freshly
# derived from the sample map present in this clone.

import argparse
import csv
from pathlib import Path

# Anchor the default paths to this script's location: both the input sample map and the output
# metadata are small TSVs that live IN the repo (no DATA_DIR needed), and resolving them here
# means the script runs correctly from any working directory.
SCRIPT_DIR = Path(__file__).resolve().parent

# Hardcoded patient-level response map (34 patients: 9 R1 / 14 R2 / 11 NR), keyed by zero-padded patient number, transcribed from GSE246613.
RESPONSE_MAP = {
    # R1 (pre-existing immunity responders)
    "02": "R1", "04": "R1", "05": "R1", "07": "R1", "09": "R1",
    "12": "R1", "16": "R1", "43": "R1", "48": "R1",
    # R2 (post-combination responders)
    "01": "R2", "03": "R2", "11": "R2", "18": "R2", "19": "R2",
    "20": "R2", "23": "R2", "30": "R2", "42": "R2", "46": "R2",
    "47": "R2", "56": "R2", "63": "R2", "66": "R2",
    # NR (non-responders)
    "06": "NR", "10": "NR", "15": "NR", "17": "NR", "26": "NR",
    "39": "NR", "45": "NR", "52": "NR", "53": "NR", "57": "NR", "64": "NR",
}

TREATMENT_MAP = {
    "A": "BASE",
    "B": "PD1",       
    "C": "RTPD1",    
}

def main():
    # Parse the input sample map path and the output filename.
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-map", default=str(SCRIPT_DIR.parent / "00_2_sample_mapping" / "sample_map_gex_final.tsv"))
    ap.add_argument("--out", default=str(SCRIPT_DIR / "sample_metadata_final.tsv"))
    args = ap.parse_args()

    # Fixed output schema (tab-separated).
    fieldnames = [
        "srr", "sample_id", "batch", "cohort",
        "treatment", "fraction", "dataset_origin", "response",
    ]

    # Counters to report how many libraries matched a response label.
    matched = 0
    unmatched = []

    with open(args.sample_map) as f_in, open(args.out, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for line in f_in:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            srr, patient, timepoint, subnum, fraction, sample_name = parts[:6]

            # Rebuild the original library label, kept as batch_key downstream (not a join key).
            batch = f"h{patient}{timepoint}{subnum}_{fraction}"

            # Join response by patient number (zero-padded defensively); flag libraries whose patient is absent from RESPONSE_MAP.
            pnum = patient.zfill(2)
            response = RESPONSE_MAP.get(pnum)
            if response is not None:
                matched += 1
            else:
                response = "NA"
                unmatched.append(batch)

            # Assemble and write one metadata row per library.
            row = {
                "srr": srr,
                "sample_id": sample_name,
                "batch": batch,
                "cohort": f"Patient{pnum}",
                "treatment": TREATMENT_MAP.get(timepoint, "NA"),
                "fraction": "CD45+" if fraction == "P" else "CD45-",
                "dataset_origin": "immune" if fraction == "P" else "non_immune",
                "response": response,
            }
            writer.writerow(row)

    # Report the match count and any unmatched libraries for validation.
    print(f"Libraries with response found: {matched}")
    if unmatched:
        print(f"Libraries without a match in RESPONSE_MAP ({len(unmatched)}): {sorted(unmatched)}")
    print(f"Output: {args.out}")


if __name__ == "__main__":
    main()
