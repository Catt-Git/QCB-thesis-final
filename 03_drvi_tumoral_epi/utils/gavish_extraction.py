"""
03 utils: split the Gavish 2023 metaprograms (GAVISH.csv, one column per MP) into one
gene list per program.

MP1 (Cell Cycle - G2/M) is missing from the CSV and from MSigDB, so it is written from the
authors' list (tiroshlab/3ca, ITH_hallmarks/MPs_distribution/MP_list.RDS.gz), with the
paper's symbols: HIST1H4C, not H4C3, which shiao.h5ad does not have.

Input : $DATA_DIR/GAVISH.csv
Output: $DATA_DIR/signatures/GAVISH_metaprograms/<program>.txt   (41 files)

Usage:
  export DATA_DIR=<data root>
  python3 gavish_extraction.py
"""
import os
from pathlib import Path

import pandas as pd

input_csv = Path(os.environ["DATA_DIR"]) / "GAVISH.csv"
output_dir = Path(os.environ["DATA_DIR"]) / "signatures" / "GAVISH_metaprograms"
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_csv)

# One gene per line, empty cells dropped.
for program in df.columns:
    genes = df[program].dropna().astype(str).str.strip()
    genes = genes[genes != ""]
    genes.to_csv(output_dir / f"{program}.txt", index=False, header=False)

MP1_CELL_CYCLE_G2_M = [
    "TOP2A", "UBE2C", "HMGB2", "NUSAP1", "CENPF", "CCNB1", "TPX2", "CKS2", "BIRC5", "PRC1",
    "PTTG1", "KPNA2", "MKI67", "CDC20", "CDK1", "CCNB2", "CDKN3", "SMC4", "NUF2", "ARL6IP1",
    "CKAP2", "ASPM", "PLK1", "CKS1B", "CCNA2", "AURKA", "MAD2L1", "GTSE1", "HMMR", "UBE2T",
    "CENPE", "CENPA", "KIF20B", "AURKB", "CDCA3", "CDCA8", "UBE2S", "KNSTRN", "KIF2C", "PBK",
    "TUBA1B", "DLGAP5", "TACC3", "STMN1", "DEPDC1", "ECT2", "CENPW", "ZWINT", "HIST1H4C",
    "KIF23",
]
(output_dir / "MP1_CELL_CYCLE_G2_M.txt").write_text("\n".join(MP1_CELL_CYCLE_G2_M) + "\n")

# The 40 CSV columns plus MP1.
expected_programs = 41
created_files = sorted(output_dir.glob("*.txt"))

print(f"programs in the CSV : {len(df.columns)}")
print(f"files written       : {len(created_files)} (expected {expected_programs})")
if len(created_files) != expected_programs:
    print("\nprograms in the CSV:")
    for program in df.columns:
        print(f"  - {program}")
    print("\nfiles written:")
    for file in created_files:
        print(f"  - {file.name}")
    raise SystemExit(f"expected {expected_programs} files")
