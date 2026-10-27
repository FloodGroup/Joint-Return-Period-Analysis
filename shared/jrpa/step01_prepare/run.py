"""01. Read observations, apply quality rules and vertical references."""

import json
import re
from pathlib import Path

from jrpa.common.analysis import Analysis
from jrpa.common.annual_samples import daily_maxima, hourly_tide_maxima
from jrpa.common.config import settings_for
from jrpa.common.datums import apply_epochs
from jrpa.common.io import (
    local_path,
    read_screened_maxima,
    read_series,
    sha256,
    write_csv,
    write_json,
)
from jrpa.common.quality import check_storm_years, exclude_intervals


def step01_prepare(manifest_path, output, settings=None):
    """01. Read sources, apply datum and quality rules, and prepare tide annual samples."""
    manifest_path = Path(manifest_path).resolve()
    root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf8"))
    settings = settings or settings_for(manifest["region"])
    if manifest["region"] != settings.region:
        raise ValueError("Manifest region differs from the requested analysis")
    if manifest.get("synthetic", False) != settings.synthetic:
        raise ValueError("Synthetic and study settings cannot be interchanged")
    out = Path(output).resolve()
    if out == root or root.is_relative_to(out):
        raise ValueError("Output directory must be separate from the private input tree")
    out.mkdir(parents=True, exist_ok=True)
    inputs = []
    tide = {}
    annual_tides = {}
    tide_quality = {}
    for gauge, item in manifest["tides"].items():
        if not re.fullmatch(r"[A-Za-z0-9_+.-]+", gauge) or gauge in {".", ".."}:
            raise ValueError("Tide identifiers must be safe path components")
        path = local_path(root, item["path"])
        s = read_series(
            path,
            "date" if settings.event_clock == "daily" else "time_utc",
            "water_m",
            daily=settings.event_clock == "daily",
        )
        inputs.append(dict(role="tide", station_id=gauge, path=item["path"], sha256=sha256(path)))
        if item.get("epochs"):
            s, audit = apply_epochs(
                s, item["epochs"], settings.datum, item.get("epoch_timezone", settings.timezone)
            )
            write_json(out / "quality" / (gauge + "_datum.json"), audit)
        elif item.get("datum") != settings.datum:
            raise ValueError(f"A documented {settings.datum} reference is required for {gauge}")
        s, audit = exclude_intervals(s, item.get("exclusions", []))
        tide[gauge] = s
        tide_quality[gauge] = audit
        annual_tides[gauge], extra = exclude_intervals(s, item.get("annual_exclusions", []))
        tide_quality[gauge].extend(extra)
    water_maxima = {
        g: (
            daily_maxima(s, "H")
            if settings.event_clock == "daily"
            else hourly_tide_maxima(s, settings.timezone, settings.annual_coverage)
        )
        for g, s in annual_tides.items()
    }
    for gauge, item in manifest["tides"].items():
        if item.get("annual_maxima"):
            path = local_path(root, item["annual_maxima"])
            water_maxima[gauge] = read_screened_maxima(path)
            inputs.append(
                dict(
                    role="screened_tide_annual_maxima",
                    station_id=gauge,
                    path=item["annual_maxima"],
                    sha256=sha256(path),
                )
            )
    if settings.region == "gba" and not settings.synthetic:
        hk = {g: s for g, s in tide.items() if g in {"TBT", "QUB", "SPW", "TPK"}}
        mainland = {g: s for g, s in tide.items() if g not in hk}
        if not hk:
            raise ValueError("GBA tide-quality analysis requires Hong Kong reference gauges")
        additions = []
        if "SZ" in mainland and not manifest.get("supplemental_peaks"):
            raise ValueError(
                "The final GBA analysis requires the private supplemental-peak table with provenance"
            )
        if manifest.get("supplemental_peaks"):
            path = local_path(root, manifest["supplemental_peaks"])
            additions = json.loads(path.read_text(encoding="utf8"))
            inputs.append(
                dict(
                    role="supplemental_observations",
                    path=manifest["supplemental_peaks"],
                    sha256=sha256(path),
                )
            )
        checked, days, response = check_storm_years(
            mainland, hk, {g: water_maxima[g] for g in mainland}, additions
        )
        water_maxima.update(checked)
        write_csv(out / "quality/regional_storm_days.csv", days)
        write_csv(out / "quality/regional_storm_responses.csv", response.reset_index())
    donors = []
    for donor in manifest.get("regional_rainfall", []):
        path = local_path(root, donor["path"])
        s = read_series(path, "date", "rain_mm", daily=True)
        maxima = daily_maxima(s, "R")
        donors.append(
            dict(
                station_id=donor["station_id"],
                location=donor["location"],
                maxima=maxima.loc[maxima.used, "maximum"].to_numpy(float),
            )
        )
        inputs.append(
            dict(
                role="regional_rainfall",
                station_id=donor["station_id"],
                path=donor["path"],
                sha256=sha256(path),
            )
        )
    return Analysis(
        manifest_path,
        manifest,
        settings,
        out,
        inputs=inputs,
        tide=tide,
        tide_quality=tide_quality,
        water_maxima=water_maxima,
        donors=donors,
    )
