"""Priority merging, rainfall event selection, and paired exposure."""

import numpy as np
import pandas as pd
from .event_kernel import event_catalog

YEAR_HOURS = 365.2425 * 24


def merge_hourly(rain_sources):
    """Prefer observed hourly values, then intact nonoverlapping 6-hour reports."""
    if not rain_sources:
        raise ValueError("No rainfall sources")
    valid_sources = [s for s in rain_sources if s["series"].notna().any()]
    if not valid_sources:
        raise ValueError("No valid rainfall observations")
    first = min(
        s["series"].index.min() - pd.Timedelta(hours=s.get("resolution_h", 1) - 1)
        for s in valid_sources
    )
    last = max(s["series"].index.max() for s in valid_sources)
    clock = pd.date_range(first, last, freq="h")
    occupied = np.zeros(len(clock), dtype=bool)
    rows = []
    ordered = sorted(valid_sources, key=lambda x: x["priority"])
    for resolution in (1, 6):
        for source in ordered:
            if source.get("resolution_h", 1) != resolution:
                continue
            values = source["series"]
            if (values.dropna() < 0).any():
                raise ValueError("Negative rainfall must be quality-controlled before merging")
            positions = clock.get_indexer(values.index)
            if (positions < 0).any():
                raise ValueError("Rainfall timestamps must lie on the hourly clock")
            if resolution == 1:
                array = values.to_numpy(float)
                take = np.isfinite(array) & ~occupied[positions]
                selected = positions[take]
                occupied[selected] = True
                rows.extend(
                    dict(
                        start_pos=int(j),
                        end_pos=int(j),
                        amount=float(v),
                        station_id=source["station_id"],
                        resolution_h=1,
                        distance_km=source.get("distance_km", 0.0),
                    )
                    for j, v in zip(selected, array[take])
                )
                continue
            for j, value in zip(positions, values.to_numpy(float)):
                left = j - resolution + 1
                if not np.isfinite(value) or left < 0 or occupied[left : j + 1].any():
                    continue
                occupied[left : j + 1] = True
                rows.append(
                    dict(
                        start_pos=left,
                        end_pos=j,
                        amount=value,
                        station_id=source["station_id"],
                        resolution_h=resolution,
                        distance_km=source.get("distance_km", 0.0),
                    )
                )
    blocks = pd.DataFrame(rows).sort_values("end_pos").reset_index(drop=True)
    return clock, blocks, occupied


def decluster(candidates, hours=72):
    if candidates.empty:
        return candidates.copy()
    groups = []
    for row in candidates.sort_values(["water_peak_utc", "selected_window_start_utc"]).to_dict(
        "records"
    ):
        if (
            not groups
            or (row["water_peak_utc"] - groups[-1][-1]["water_peak_utc"]).total_seconds() / 3600
            >= hours
        ):
            groups.append([])
        groups[-1].append(row)
    selected = []
    for group in groups:
        row = max(
            group,
            key=lambda r: (
                r["selected_rain_mm"],
                -pd.Timestamp(r["selected_window_start_utc"]).value,
            ),
        )
        selected.append(
            {
                **row,
                "cluster_members": len(group),
                "rain_mm": row["selected_rain_mm"],
                "water_max_m": row["paired_water_max_m"],
            }
        )
    return pd.DataFrame(selected)


def hourly_events(unit, rain_sources, tide, settings):
    clock, blocks, observed = merge_hourly(rain_sources)
    overlap = (clock >= tide.index.min()) & (clock <= tide.index.max())
    if not overlap.any():
        raise ValueError("Rainfall and tide records do not overlap")
    start = int(np.flatnonzero(overlap)[0])
    end = int(np.flatnonzero(overlap)[-1])
    blocks = blocks.loc[(blocks.start_pos >= start) & (blocks.end_pos <= end)].copy()
    blocks[["start_pos", "end_pos"]] -= start
    clock = clock[start : end + 1]
    observed = np.zeros(len(clock), bool)
    for r in blocks.itertuples():
        observed[r.start_pos : r.end_pos + 1] = True
    water = tide.reindex(clock).to_numpy(float)
    all_events = event_catalog(unit, blocks, clock, water, [])
    if all_events.empty or not all_events.selected_rain_mm.notna().any():
        raise ValueError("No complete rainfall event windows")
    # Threshold before water-coverage rejection and temporal declustering.
    r95 = float(all_events.selected_rain_mm.dropna().quantile(settings.rainfall_quantile))
    candidates = all_events.loc[
        (all_events.selected_rain_mm > r95) & all_events.pair_status.eq("valid")
    ].copy()
    pot = decluster(candidates, settings.decluster_hours)
    exposure = observed & np.isfinite(water)
    years = float(exposure.sum() / YEAR_HOURS)
    if years <= 0:
        raise ValueError("No simultaneous observed exposure")
    meta = dict(
        unit=unit,
        R_threshold_mm=r95,
        effective_years=years,
        events_per_year=len(pot) / years,
        paired_observed_hours=int(exposure.sum()),
    )
    return all_events, pot, blocks, meta


def daily_events(unit, rain, tide, settings, excluded_dates=(), analysis_span=None):
    if rain.dropna().empty or tide.dropna().empty:
        raise ValueError("Finite rainfall and water-level observations required")
    # Retain missing endpoint dates in a declared analysis clock: they still
    # belong to the threshold's rainfall population, but never to exposure.
    start = max(rain.index.min(), tide.index.min())
    end = min(rain.index.max(), tide.index.max())
    if analysis_span is not None:
        start, end = map(pd.Timestamp, analysis_span)
    frame = pd.concat([rain.rename("rain_mm"), tide.rename("water_max_m")], axis=1).reindex(
        pd.date_range(start, end, freq="D")
    )
    if (frame.rain_mm.dropna() < 0).any():
        raise ValueError("Negative rainfall is invalid")
    positive = frame.loc[frame.rain_mm.gt(0), "rain_mm"]
    if positive.empty:
        raise ValueError("No positive rainfall")
    r95 = float(positive.quantile(settings.rainfall_quantile))
    valid = frame.notna().all(axis=1)
    candidates = frame.loc[valid & frame.rain_mm.gt(r95)].copy()
    groups = []
    for date, row in candidates.iterrows():
        if not groups or (date - groups[-1][-1][0]).days >= 3:
            groups.append([])
        groups[-1].append((date, row))
    rows = []
    for i, group in enumerate(groups):
        date, row = max(group, key=lambda x: (x[1].rain_mm, x[1].water_max_m))
        if date.normalize() in pd.to_datetime(list(excluded_dates)):
            continue
        rows.append(
            dict(
                event_id=f"{unit}_{i + 1:04d}",
                date=date,
                rain_mm=round(float(row.rain_mm), 3),
                water_max_m=round(float(row.water_max_m), 3),
                cluster_members=len(group),
            )
        )
    years = float(valid.sum() / 365.2425)
    if years <= 0:
        raise ValueError("No simultaneous daily observations")
    return (
        frame,
        pd.DataFrame(rows),
        dict(
            unit=unit,
            R_threshold_mm=r95,
            effective_years=years,
            events_per_year=len(rows) / years,
            paired_observed_days=int(valid.sum()),
        ),
    )


def merge_daily(sources):
    ordered = sorted(sources, key=lambda s: s["priority"])
    index = ordered[0]["series"].index
    for source in ordered[1:]:
        index = index.union(source["series"].index)
    frame = pd.DataFrame(index=index.sort_values(), columns=["rain_mm", "station_id"])
    for source in ordered:
        s = source["series"].reindex(frame.index)
        take = frame.rain_mm.isna() & s.notna()
        frame.loc[take, "rain_mm"] = s.loc[take]
        frame.loc[take, "station_id"] = source["station_id"]
    frame["rain_mm"] = pd.to_numeric(frame.rain_mm)
    return frame
