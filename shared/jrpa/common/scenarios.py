"""Ninety design cases; numerical climate inputs are supplied privately."""

import pandas as pd
from .config import settings_for

CLIMATES = ("C", "S126_2060", "S126_2100", "S245_2060", "S245_2100", "S370_2060", "S370_2100")
PERIODS = (10, 50, 100, 200)


def case_definitions(region):
    peak = settings_for(region).rain_peak_hour
    cases = []
    for hi, hazard in enumerate(("R", "S", "RS")):
        for ci, climate in enumerate(CLIMATES):
            for ti, period in enumerate(PERIODS):
                number = hi * 28 + ci * 4 + ti + 1
                code = {
                    "R": f"R{period:03d}",
                    "S": f"S{period:03d}",
                    "RS": f"R{period:03d}S{period:03d}",
                }[hazard]
                cases.append(
                    dict(
                        case=f"C{number:02d}",
                        name=f"C{number:02d}_{code}_{climate}",
                        hazard=hazard,
                        climate=climate,
                        return_period_years=period,
                        rain_peak_hour=peak if hazard != "S" else None,
                        tide_peak_hour=peak if hazard != "R" else None,
                        base_case=None,
                    )
                )
    for i, offset in enumerate((-9, -6, -3, 3, 6, 9), 85):
        cases.append(
            dict(
                case=f"C{i:02d}",
                name=f"C{i:02d}_R100S100_S370_2100_{int(offset / 12 * 100)}",
                hazard="RS",
                climate="S370_2100",
                return_period_years=100,
                rain_peak_hour=peak,
                tide_peak_hour=peak + offset,
                base_case="C83",
            )
        )
    return cases


def case_designs(region, designs, unit_gauges, climate_inputs):
    """Apply rainfall multipliers and model-datum sea-level increments."""
    rows = []
    d = pd.DataFrame(designs)
    for case in case_definitions(region):
        climate = climate_inputs[case["climate"]]
        for unit, gauge in unit_gauges.items():
            hazard = case["hazard"]
            selected = d.loc[d.return_period_years.eq(case["return_period_years"])]
            if hazard == "S":
                selected = selected.loc[selected.scenario.eq("S2") & selected.gauge.eq(gauge)]
            else:
                scenario = "S1" if hazard == "R" else "Kendall"
                selected = selected.loc[selected.scenario.eq(scenario) & selected.unit.eq(unit)]
            if len(selected) != 1:
                raise ValueError(f"One design record required for {case['case']}/{unit}")
            point = selected.iloc[0]
            rainfall = (
                0.0 if hazard == "S" else float(point.R_mm) * float(climate["rainfall_factor"])
            )
            water = None if hazard == "R" else float(point.H_model_m) + float(climate["SLR_m"])
            rows.append(
                dict(
                    **case,
                    unit=unit,
                    gauge=gauge,
                    R_mm=rainfall,
                    H_peak_model_m=water,
                    diagnostic_only=bool(point.diagnostic_only),
                )
            )
    return pd.DataFrame(rows)
