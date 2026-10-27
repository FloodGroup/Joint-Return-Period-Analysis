"""Vertical-reference conversion with explicit station and epoch provenance."""

import numpy as np
from jrpa.gba.settings import GBA_MAINLAND_MODEL_OFFSET, GBA_HK_MODEL_OFFSET
import pandas as pd


def apply_epochs(series, epochs, target_datum, timezone="UTC"):
    """Apply additive epoch offsets; unknown epochs remain missing."""
    if not isinstance(series.index, pd.DatetimeIndex) or series.index.tz is None:
        raise ValueError("Datum conversion requires timezone-aware observations")
    result = pd.Series(np.nan, index=series.index, dtype=float)
    assigned = pd.Series(False, index=series.index)
    audit = []
    for epoch in epochs:
        if epoch["target_datum"] != target_datum or not epoch.get("source"):
            raise ValueError("Each datum epoch needs its target datum and evidence source")
        start = pd.Timestamp(epoch["start"])
        end = pd.Timestamp(epoch["end"])
        start = start.tz_localize(timezone) if start.tzinfo is None else start
        end = end.tz_localize(timezone) if end.tzinfo is None else end
        mask = (series.index >= start) & (series.index < end)
        if (assigned & mask).any():
            raise ValueError("Overlapping datum epochs")
        assigned |= mask
        offset = epoch.get("offset_m")
        if offset is not None:
            result.loc[mask] = series.loc[mask] + float(offset)
        audit.append({**epoch, "observations": int(series.loc[mask].notna().sum())})
    return result, audit


def model_offset(region, gauge):
    if region != "gba":
        return 0.0
    return (
        GBA_HK_MODEL_OFFSET if gauge in {"TBT", "QUB", "SPW", "TPK"} else GBA_MAINLAND_MODEL_OFFSET
    )
