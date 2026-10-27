"""Deterministic artificial observations for executable public examples."""

from pathlib import Path
import numpy as np
import pandas as pd
from .config import settings_for
from .io import write_csv, write_json


def create_demo(directory, region="tokyo", years=16, seed=4217):
    if years < 8:
        raise ValueError("The demonstration needs at least eight years")
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=True)
    settings = settings_for(region)
    rng = np.random.default_rng(seed)
    if settings.event_clock == "daily":
        clock = pd.date_range("1990-01-01", periods=years * 365, freq="D")
        rain = np.where(rng.random(len(clock)) < 0.35, rng.gamma(1.4, 19, len(clock)), 0.0)
        water = 0.8 + 0.08 * np.log1p(rain) + rng.normal(0, 0.18, len(clock))
        time = "date"
    else:
        clock = pd.date_range("1990-10-01", periods=years * 365 * 24, freq="h", tz="UTC")
        rain = np.zeros(len(clock))
        starts = np.flatnonzero(rng.random(len(clock)) < 0.014)
        for start in starts:
            length = int(rng.integers(2, 20))
            end = min(start + length, len(clock))
            rain[start:end] += rng.gamma(1.5, 2.8, end - start)
        water = (
            0.8
            + 0.55 * np.sin(np.arange(len(clock)) * 2 * np.pi / 12.42)
            + 0.07 * np.log1p(rain)
            + rng.normal(0, 0.12, len(clock))
        )
        time = "time_utc"
    alternate = rain.copy()
    primary = rain.copy()
    primary[rng.random(len(clock)) < 0.025] = np.nan
    write_csv(root / "rain_primary.csv", pd.DataFrame({time: clock, "rain_mm": primary}))
    write_csv(root / "rain_alternate.csv", pd.DataFrame({time: clock, "rain_mm": alternate}))
    write_csv(root / "tide.csv", pd.DataFrame({time: clock, "water_m": water}))
    manifest = dict(
        region=region,
        synthetic=True,
        tides={"SYN_T01": dict(path="tide.csv", datum=settings.datum)},
        units=[
            dict(
                unit="SYN01",
                gauge="SYN_T01",
                rainfall_sources=[
                    dict(
                        station_id="SYN_R01",
                        path="rain_primary.csv",
                        priority=1,
                        resolution_h=24 if settings.event_clock == "daily" else 1,
                    ),
                    dict(
                        station_id="SYN_R02",
                        path="rain_alternate.csv",
                        priority=2,
                        resolution_h=24 if settings.event_clock == "daily" else 1,
                    ),
                ],
            )
        ],
    )
    write_json(root / "manifest.json", manifest)
    return root / "manifest.json"
