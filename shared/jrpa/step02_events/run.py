"""02. Merge rainfall sources and select independent event pairs."""

import re

from jrpa.common.annual_samples import daily_maxima, hourly_rain_maxima
from jrpa.common.io import (
    local_path,
    read_screened_maxima,
    read_series,
    sha256,
    write_csv,
    write_json,
)
from jrpa.common.quality import exclude_intervals
from jrpa.gba.settings import GBA_EVENT_EXCLUSIONS

from .coverage import coverage_figure
from .events import daily_events, hourly_events, merge_daily, merge_hourly


def step02_events(analysis):
    """02. Merge rainfall, build paired events, and prepare rainfall annual samples."""
    analysis.require(2)
    settings, manifest, out = analysis.settings, analysis.manifest, analysis.out
    root = analysis.manifest_path.parent
    inputs, tide = analysis.inputs, analysis.tide
    sources = []
    seen = set()
    for unit in manifest["units"]:
        uid = unit["unit"]
        gauge = unit["gauge"]
        if not re.fullmatch(r"[A-Za-z0-9_-]+", uid) or uid in seen:
            raise ValueError("Unit identifiers must be unique safe path components")
        seen.add(uid)
        loaded = []
        annual_loaded = []
        for source in unit["rainfall_sources"]:
            if not re.fullmatch(r"[A-Za-z0-9_-]+", source["station_id"]):
                raise ValueError("Rainfall identifiers must be safe path components")
            path = local_path(root, source["path"])
            s = read_series(
                path,
                "date" if settings.event_clock == "daily" else "time_utc",
                "rain_mm",
                daily=settings.event_clock == "daily",
            )
            s, audit = exclude_intervals(s, source.get("exclusions", []))
            loaded.append({**source, "series": s})
            annual_s, _ = exclude_intervals(s, source.get("annual_exclusions", []))
            annual_loaded.append({**source, "series": annual_s})
            inputs.append(
                dict(
                    role="rainfall",
                    station_id=source["station_id"],
                    path=source["path"],
                    sha256=sha256(path),
                )
            )
            write_json(out / "units" / uid / (source["station_id"] + "_quality.json"), audit)
        folder = out / "units" / uid
        if settings.event_clock == "daily":
            merged = merge_daily(loaded)
            coverage_figure(folder, merged.index, daily_sources=merged.station_id)
            all_events, pot, event_meta = daily_events(
                uid,
                merged.rain_mm,
                tide[gauge],
                settings,
                GBA_EVENT_EXCLUSIONS.get(uid, []) if not settings.synthetic else [],
                analysis_span=unit.get("event_span"),
            )
            rain_maxima = daily_maxima(merged.rain_mm, "R")
            for sid, count in (
                merged.loc[merged.rain_mm.notna(), "station_id"].value_counts().items()
            ):
                sources.append(dict(unit=uid, station_id=sid, selected_observations=int(count)))
        else:
            all_events, pot, blocks, event_meta = hourly_events(uid, loaded, tide[gauge], settings)
            fullclock, fullblocks, _ = merge_hourly(annual_loaded)
            coverage_figure(folder, fullclock, blocks=fullblocks)
            rain_maxima = hourly_rain_maxima(
                fullclock, fullblocks, settings.timezone, settings.annual_coverage
            )
            write_csv(folder / "rainfall_intervals.csv", fullblocks)
            for sid, group in fullblocks.groupby("station_id"):
                sources.append(
                    dict(
                        unit=uid,
                        station_id=sid,
                        selected_observations=int(group.resolution_h.sum()),
                    )
                )
        if unit.get("annual_maxima"):
            path = local_path(root, unit["annual_maxima"])
            rain_maxima = read_screened_maxima(path)
            inputs.append(
                dict(
                    role="screened_rain_annual_maxima",
                    unit=uid,
                    path=unit["annual_maxima"],
                    sha256=sha256(path),
                )
            )
        event_meta.update(gauge=gauge, datum=settings.datum)
        write_csv(folder / "all_events.csv", all_events.reset_index())
        write_csv(folder / "pot_events.csv", pot)
        write_json(folder / "event_settings.json", event_meta)
        write_csv(folder / "rain_annual_maxima.csv", rain_maxima)
        analysis.units.append(
            dict(unit=unit, folder=folder, pot=pot, event_meta=event_meta, rain_maxima=rain_maxima)
        )
    analysis.sources = sources
    analysis.completed = 2
