#!/usr/bin/env python3
"""Print `<run_id> <compartment> <out_root>` for the current environment, one line.

Used by `signature_interpretation_all.sh`, so the run id (`cell_set`) and the output root
(`signature_common`, $OUT_ROOT) are not re-derived in bash.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import cell_set as CS              # noqa: E402
import signature_common as C       # noqa: E402

print(C.RUN_ID, CS.compartment(), C.OUT_ROOT)
