"""Strict local input and portable result serialization."""

import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean_json(v) for v in value]
    if isinstance(value, np.generic):
        return clean_json(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, (pd.Timestamp, Path)):
        return str(value)
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(clean_json(value), indent=2, allow_nan=False) + "\n", encoding="utf8"
    )


def write_csv(path, frame):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(frame).to_csv(path, index=False, encoding="utf8", float_format="%.15g")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def read_series(path, time_column, value_column, daily=False):
    frame = pd.read_csv(path)
    if not {time_column, value_column} <= set(frame):
        raise ValueError(f"Required columns: {time_column}, {value_column}")
    index = pd.to_datetime(frame[time_column], utc=not daily, errors="raise")
    if index.duplicated().any():
        raise ValueError(f"Duplicate timestamps: {Path(path).name}")
    values = pd.to_numeric(frame[value_column], errors="raise")
    if np.isinf(values.to_numpy(float)).any():
        raise ValueError("Infinite observations are invalid")
    timestamps = pd.DatetimeIndex(index)
    if (
        not daily
        and (
            (timestamps.minute != 0) | (timestamps.second != 0) | (timestamps.microsecond != 0)
        ).any()
    ):
        raise ValueError("Hourly observations must use exact UTC clock hours")
    if daily and not timestamps.equals(timestamps.normalize()):
        raise ValueError("Daily observations must use date labels at midnight")
    return pd.Series(values.to_numpy(float), index=index).sort_index()


def local_path(root, relative):
    """Manifest paths must remain within the explicitly supplied private root."""
    root = Path(root).resolve()
    target = (root / relative).resolve()
    if not target.is_relative_to(root):
        raise ValueError("Manifest path escapes its private input directory")
    return target


def read_screened_maxima(path):
    """Read privately archived annual observations, never fitted design levels."""
    frame = pd.read_csv(path)
    if not {"year", "maximum", "used"} <= set(frame):
        raise ValueError("Screened annual samples require year, maximum, and used")
    if frame.year.duplicated().any() or not frame.used.isin([True, False]).all():
        raise ValueError("Unique years and explicit Boolean selection flags required")
    if not np.isfinite(frame.loc[frame.used, "maximum"]).all():
        raise ValueError("Accepted annual maxima must be finite")
    return frame.sort_values("year").reset_index(drop=True)
