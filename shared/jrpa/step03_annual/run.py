"""03. Fit annual rainfall (S1) and water-level (S2) distributions."""

import numpy as np

from jrpa.common.config import tide_transfers
from jrpa.common.datums import model_offset
from jrpa.common.io import write_csv, write_json

from .annual import fit_annual, transfer_curve
from .regional import regional_design


def step03_annual(analysis):
    """03. Fit annual S1 rainfall and S2 water-level frequency curves."""
    analysis.require(3)
    settings, out = analysis.settings, analysis.out
    water_maxima, tide, donors = analysis.water_maxima, analysis.tide, analysis.donors
    annual_tide = {}
    tide_design = {}
    design_rows = []
    for gauge, maxima in water_maxima.items():
        # All annual tide fits use the full station record, independently of rainfall.
        meta, values, candidates = fit_annual(maxima, "H", settings)
        annual_tide[gauge] = meta
        tide_design[gauge] = values
        write_csv(out / "tides" / gauge / "annual_maxima.csv", maxima)
        write_csv(out / "tides" / gauge / "candidates.csv", candidates)
    for short, long in tide_transfers(settings.region).items():
        if short in tide:
            if long not in tide or not tide_design[long]:
                raise ValueError(f"Transfer donor required: {long}")
            tide_design[short], audit = transfer_curve(
                water_maxima[short],
                water_maxima[long],
                tide_design[long],
                minimum_year=2019 if settings.region == "san_francisco" else None,
            )
            annual_tide[short]["transfer"] = dict(donor=long, **audit)
    for gauge, values in tide_design.items():
        for row in values:
            design_rows.append(
                dict(
                    gauge=gauge,
                    scenario="S2",
                    return_period_years=row["return_period_years"],
                    R_mm=0.0,
                    H_m=row["level"],
                    H_model_m=row["level"] + model_offset(settings.region, gauge),
                    diagnostic_only=row["diagnostic_only"],
                    synthetic=settings.synthetic,
                )
            )
    for item in analysis.units:
        unit, folder, rain_maxima = item["unit"], item["folder"], item["rain_maxima"]
        uid, gauge = unit["unit"], unit["gauge"]
        if (
            settings.region == "gba"
            and uid in {"01", "02", "03a", "03b", "04", "05", "06", "07", "08"}
            and not settings.synthetic
        ):
            values, annual_meta = regional_design(
                rain_maxima.loc[rain_maxima.used, "maximum"].to_numpy(float),
                unit["rain_location"],
                donors,
                settings.periods,
                settings.rainfall_factor_24h,
            )
            annual_meta["diagnostic_only"] = False
            for row in values:
                row["diagnostic_only"] = False
            candidates = []
        else:
            annual_meta, values, candidates = fit_annual(rain_maxima, "R", settings)
            for row in values:
                row["level"] *= settings.rainfall_factor_24h
        write_json(folder / "rain_annual_model.json", annual_meta)
        write_csv(folder / "rain_annual_candidates.csv", candidates)
        for row in values:
            design_rows.append(
                dict(
                    unit=uid,
                    gauge=gauge,
                    scenario="S1",
                    return_period_years=row["return_period_years"],
                    R_mm=row["level"],
                    H_m=np.nan,
                    H_model_m=np.nan,
                    diagnostic_only=row["diagnostic_only"],
                    synthetic=settings.synthetic,
                )
            )
    analysis.annual_tide = annual_tide
    analysis.design_rows = design_rows
    analysis.completed = 3
