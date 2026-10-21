"""
02 utils - compatibility of scib 1.1.7 with numpy 2.x / pandas 3.x

Import BEFORE scib to restore APIs removed in numpy 2.4 / pandas 3.0, which
otherwise crash scib hours into a metrics job:

    import scib_compat  # noqa: F401  (must precede the scib import)
    import scib

Restores np.in1d, pd.value_counts and positional Series[int] access, and fixes
anndata2ri.activate/deactivate. Each patch applies only if the API is missing.
"""

import os
import shutil

import numpy as np
import pandas as pd

# numpy 

if not hasattr(np, "in1d"):
    np.in1d = np.isin

# pandas

if not hasattr(pd, "value_counts"):

    def _value_counts(values, *args, **kwargs):
        """Equivalent of the module-level function removed in pandas 3.0."""
        series = values if isinstance(values, pd.Series) else pd.Series(values)
        return series.value_counts(*args, **kwargs)

    pd.value_counts = _value_counts

# rpy2 / anndata2ri


def require_r():
    """Stop with a readable message if R is unreachable (env not activated)."""
    if not os.environ.get("R_HOME") and shutil.which("R") is None:
        raise SystemExit(
            "R is not reachable: rpy2 cannot start, so anndata2ri cannot be "
            "imported and this module fails before any metric runs.\n"
            "Activate the environment before launching the script:\n"
            "    conda activate benchmark-py-r   # locally and on the cluster"
        )


require_r()

import anndata2ri  # noqa: E402
from rpy2.robjects import conversion, default_converter  # noqa: E402


def _activate():
    conversion.set_conversion(default_converter + anndata2ri.converter)


def _deactivate():
    conversion.set_conversion(default_converter)


anndata2ri.activate = _activate
anndata2ri.deactivate = _deactivate

# pandas: positional access on Series
#
# scib uses series[0] on non-integer indexes (KeyError in pandas 3); fall back to
# iloc only in that case.

_series_getitem = pd.Series.__getitem__


def _series_getitem_positional_fallback(self, key):
    try:
        return _series_getitem(self, key)
    except KeyError:
        if isinstance(key, (int, np.integer)) and not pd.api.types.is_integer_dtype(
            self.index
        ):
            return self.iloc[key]
        raise


pd.Series.__getitem__ = _series_getitem_positional_fallback
