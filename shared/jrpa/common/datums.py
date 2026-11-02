"""Vertical-reference conversion with explicit station and epoch provenance."""

import numpy as np
from jrpa.gba.settings import GBA_MAINLAND_MODEL_OFFSET, GBA_HK_MODEL_OFFSET
import pandas as pd


def apply_epochs(series, epochs, target_datum, timezone="UTC", *, daily=False):
    """Apply additive epoch offsets; unknown epochs remain missing."""
    if not isinstance(series.index, pd.DatetimeIndex) or series.index.hasnans:
        raise ValueError("Datum conversion requires valid observation dates")
    if daily:
        if series.index.tz is not None or not series.index.equals(series.index.normalize()):
            raise ValueError("Daily observations require local date labels without a timezone")
    elif series.index.tz is None:
        raise ValueError("Datum conversion requires timezone-aware observations")
    result = pd.Series(np.nan, index=series.index, dtype=float)
    assigned = pd.Series(False, index=series.index)
    audit = []
    intervals = []
    for epoch in epochs:
        if epoch["target_datum"] != target_datum or not epoch.get("source"):
            raise ValueError("Each datum epoch needs its target datum and evidence source")
        start = pd.Timestamp(epoch["start"])
        end = pd.Timestamp(epoch["end"])
        start = start.tz_localize(timezone) if start.tzinfo is None else start
        end = end.tz_localize(timezone) if end.tzinfo is None else end
        if pd.isna(start) or pd.isna(end) or start >= end:
            raise ValueError("Datum epochs require a valid start before the exclusive end")
        if daily:
            start = start.tz_convert(timezone).tz_localize(None)
            end = end.tz_convert(timezone).tz_localize(None)
            if start != start.normalize() or end != end.normalize():
                raise ValueError("Daily datum epochs must start and end at local midnight")
        if any(
            start < previous_end and end > previous_start
            for previous_start, previous_end in intervals
        ):
            raise ValueError("Overlapping datum epochs")
        intervals.append((start, end))
        mask = (series.index >= start) & (series.index < end)
        if (assigned & mask).any():
            raise ValueError("Overlapping datum epochs")
        assigned |= mask
        offset = epoch.get("offset_m")
        if offset is not None:
            if not np.isfinite(float(offset)):
                raise ValueError("Datum offsets must be finite or null for unknown offsets")
            result.loc[mask] = series.loc[mask] + float(offset)
        audit.append({**epoch, "observations": int(series.loc[mask].notna().sum())})
    return result, audit


def model_offset(region, gauge):
    if region != "gba":
        return 0.0
    return (
        GBA_HK_MODEL_OFFSET if gauge in {"TBT", "QUB", "SPW", "TPK"} else GBA_MAINLAND_MODEL_OFFSET
    )
