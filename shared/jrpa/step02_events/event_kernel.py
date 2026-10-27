"""Accumulation-interval event pairing without rainfall disaggregation."""

import math
import numpy as np
import pandas as pd

DRY_H = 12


def event_catalog(unit, blocks, clock, h, inventory_rows):
    n = len(clock)
    coverage = np.zeros(n, bool)
    known_dry = np.zeros(n, bool)
    end_amount = np.zeros(n)
    amounts = blocks.amount.to_numpy(float)
    left = blocks.start_pos.to_numpy(int)
    right = blocks.end_pos.to_numpy(int)
    for a, b, v in zip(left, right, amounts):
        coverage[a : b + 1] = True
        if v == 0:
            known_dry[a : b + 1] = True
        end_amount[b] = v
    cumamount = np.r_[0, np.cumsum(end_amount)]
    cumcoverage = np.r_[0, np.cumsum(coverage)]
    wet = np.flatnonzero(amounts > 0)
    if not len(wet):
        raise ValueError("No wet accumulation intervals")
    separate = np.r_[True, left[wet[1:]] - right[wet[:-1]] - 1 >= DRY_H]
    starts = wet[separate]
    lasts = wet[np.r_[np.flatnonzero(separate)[1:] - 1, len(wet) - 1]]
    dry12 = pd.Series(known_dry).rolling(12, min_periods=12).sum().eq(12).to_numpy()
    rows = []
    for k, (a, b) in enumerate(zip(starts, lasts), 1):
        start, end = (int(left[a]), int(right[b]))
        span = end - start + 1
        prior = int(right[lasts[k - 2]]) if k > 1 else None
        dryconfirmed = bool(dry12[prior + 12 : start].any()) if prior is not None else None
        rec = dict(
            event_id=f"{unit[:3]}_rain_{k:05d}",
            event_start_utc=clock[start] - pd.Timedelta(hours=1),
            event_end_utc=clock[end],
            event_duration_h=span,
            prior_dry_gap_confirmed=dryconfirmed,
            selected_rain_mm=np.nan,
            pair_status="no_complete_rain_window",
        )
        if span < 24:
            window_start, window_end = (start, end)
            if cumcoverage[end + 1] - cumcoverage[start] != span:
                rows.append(rec)
                continue
        else:
            ends = right[(right >= start + 23) & (right <= end)]
            starts24 = ends - 23
            whole_left = np.isin(starts24, left)
            full = cumcoverage[ends + 1] - cumcoverage[starts24] == 24
            ends = ends[whole_left & full]
            starts24 = starts24[whole_left & full]
            if not len(ends):
                rows.append(rec)
                continue
            totals = cumamount[ends + 1] - cumamount[starts24]
            ix = int(np.argmax(totals))
            window_start = int(starts24[ix])
            window_end = int(ends[ix])
        rain = float(cumamount[window_end + 1] - cumamount[window_start])
        hh = h[window_start : window_end + 1]
        nh = int(np.isfinite(hh).sum())
        nw = len(hh)
        part = blocks[(blocks.start_pos >= window_start) & (blocks.end_pos <= window_end)]
        source_ids = "|".join(sorted(set(part.station_id)))
        watermax = float(np.nanmax(hh)) if nh else np.nan
        peak = clock[window_start + int(np.nanargmax(hh))] if nh else pd.NaT
        rec.update(
            selected_window_start_utc=clock[window_start] - pd.Timedelta(hours=1),
            selected_window_end_utc=clock[window_end],
            selected_window_h=nw,
            selected_rain_mm=rain,
            paired_water_max_m=watermax,
            water_peak_utc=peak,
            water_valid_hours=nh,
            water_coverage_fraction=nh / nw,
            rain_source_station_ids=source_ids,
            rain_source_count=len(set(part.station_id)),
            uses_6h=bool(part.resolution_h.eq(6).any()),
            maximum_source_distance_km=float(part.distance_km.max()),
            selected_rule="whole_event" if span < 24 else "max_24h",
            pair_status="valid" if nh >= math.ceil(0.75 * nw) else "water_coverage_below_75pct",
        )
        rows.append(rec)
    return pd.DataFrame(rows)
