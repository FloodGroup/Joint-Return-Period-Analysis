"""Annual-maximum samples. Fitting these samples is step 03."""

import numpy as np
import pandas as pd


def hourly_rain_maxima(clock, blocks, timezone, coverage_gate=0.75, boundary_deduplicate=True):
    occupied = np.zeros(len(clock), bool)
    totals = np.zeros(len(clock))
    for row in blocks.itertuples():
        occupied[row.start_pos : row.end_pos + 1] = True
        totals[row.end_pos] = row.amount
    cs = np.r_[0, np.cumsum(totals)]
    cv = np.r_[0, np.cumsum(occupied)]
    ends = blocks.end_pos.to_numpy(int)
    ends = ends[ends >= 23]
    starts = ends - 23
    valid = np.isin(starts, blocks.start_pos.to_numpy(int)) & (cv[ends + 1] - cv[starts] == 24)
    ends = ends[valid]
    starts = starts[valid]
    amount = cs[ends + 1] - cs[starts]
    local = clock.tz_convert(timezone)
    year = local.year + (local.month >= 10)
    if boundary_deduplicate:
        chosen = {}
        selected_windows = []
        for i in np.argsort(-amount, kind="stable"):
            wy = int(year[ends[i]])
            if wy in chosen:
                continue
            left, right = int(starts[i]), int(ends[i])
            if any(left <= b and right >= a for a, b in selected_windows):
                continue
            chosen[wy] = float(amount[i])
            selected_windows.append((left, right))
            if len(chosen) == len(np.unique(year[ends])):
                break
        values = pd.Series(chosen, dtype=float)
    else:
        values = pd.Series(amount, index=year[ends]).groupby(level=0).max()
    counts = pd.Series(occupied, index=year).groupby(level=0).sum()
    rows = []
    for wy, count in counts.items():
        expected = (
            pd.Timestamp(int(wy), 10, 1) - pd.Timestamp(int(wy) - 1, 10, 1)
        ).total_seconds() / 3600
        coverage = count / expected
        rows.append(
            dict(
                year=int(wy),
                maximum=values.get(wy, np.nan),
                coverage=coverage,
                used=bool(coverage >= coverage_gate and np.isfinite(values.get(wy, np.nan))),
            )
        )
    return pd.DataFrame(rows)

def hourly_tide_maxima(series, timezone, coverage_gate=0.75):
    local = series.index.tz_convert(timezone)
    year = local.year + (local.month >= 10)
    frame = pd.DataFrame({"year": year, "value": series.to_numpy(float)})
    rows = []
    for wy, part in frame.groupby("year"):
        expected = (
            pd.Timestamp(int(wy), 10, 1) - pd.Timestamp(int(wy) - 1, 10, 1)
        ).total_seconds() / 3600
        coverage = part.value.notna().sum() / expected
        rows.append(
            dict(
                year=int(wy),
                maximum=float(part.value.max()),
                coverage=coverage,
                used=bool(coverage >= coverage_gate and part.value.notna().any()),
            )
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["complete"] = result.used
        median = result.loc[result.complete, "maximum"].median()
        result["partial_included"] = (~result.complete) & result.maximum.ge(median)
        result["used"] = result.complete | result.partial_included
    return result

def daily_maxima(series, variable):
    """Calendar years; GBA tide coverage retains the declared seasonal/partial rule."""
    s = series.dropna()
    rows = []
    for year, part in s.groupby(s.index.year):
        coverage = len(part) / 365.25
        jas = part[(part.index.month >= 7) & (part.index.month <= 9)]
        rows.append(
            dict(
                year=int(year),
                maximum=float(part.max()),
                coverage=coverage,
                coverage_JAS=len(jas) / 92,
                used=coverage >= 0.75 if variable == "H" else len(part) >= 274,
            )
        )
    frame = pd.DataFrame(rows)
    if variable == "H" and not frame.empty:
        complete = (frame.coverage >= 0.75) | (frame.coverage_JAS >= 0.8)
        median = frame.loc[complete, "maximum"].median()
        frame["partial_included"] = (~complete) & (frame.maximum >= median)
        frame["used"] = complete | frame.partial_included
    return frame
