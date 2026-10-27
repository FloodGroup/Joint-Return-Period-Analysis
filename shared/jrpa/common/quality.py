"""Explicit exclusions and storm-day checks, with observation values kept private."""

import numpy as np
import pandas as pd


def exclude_intervals(series, exclusions):
    result = series.copy()
    audit = []
    for record in exclusions:
        if not record.get("reason"):
            raise ValueError("Every exclusion needs a reason")
        start = pd.Timestamp(record["start"])
        end = pd.Timestamp(record.get("end", record["start"]))
        if result.index.tz is not None:
            start = start.tz_localize("UTC") if start.tzinfo is None else start
            end = end.tz_localize("UTC") if end.tzinfo is None else end
        mask = (result.index >= start) & (result.index <= end)
        audit.append({**record, "removed": int(result.loc[mask].notna().sum())})
        result.loc[mask] = np.nan
    return result, audit


def storm_days(hong_kong):
    anomalies = []
    for station, series in hong_kong.items():
        s = series.dropna()
        anomalies.append((s - s.groupby(s.index.year).transform("median")).rename(station))
    combined = pd.concat(anomalies, axis=1)
    peak = combined.max(axis=1)
    retained = []
    for date, value in peak[peak >= 1.0].sort_values(ascending=False).items():
        if all(abs((date - other).days) > 3 for other, _ in retained):
            retained.append((date, value))
    return pd.DataFrame(retained, columns=["date", "HK_anomaly_max_m"]).sort_values("date")


def check_storm_years(mainland, hong_kong, maxima, observed_additions=()):
    days = storm_days(hong_kong)
    responses = {}
    for station, series in mainland.items():
        s = series.dropna()
        median = s.groupby(s.index.year).median()
        responses[station] = {
            row.date: float(
                s.loc[row.date - pd.Timedelta(days=1) : row.date + pd.Timedelta(days=1)].max()
            )
            - float(median.get(row.date.year, np.nan))
            for row in days.itertuples()
        }
    responses = pd.DataFrame(responses)
    result = {k: v.copy() for k, v in maxima.items()}
    for station, table in result.items():
        table["storm_quality_flag"] = ""
        for row in days.itertuples():
            others = responses.loc[row.date].drop(station).dropna()
            est = float(others.median()) if len(others) else np.nan
            strong = est >= 0.8 if len(others) else row.HK_anomaly_max_m >= 1.5
            own = responses.loc[row.date, station]
            if strong and (not np.isfinite(own) or own < est - 0.8):
                mask = table.year.eq(row.date.year)
                table.loc[mask, "used"] = False
                table.loc[mask, "storm_quality_flag"] = (
                    "Missing or distorted regional storm-tide peak"
                )
    # These are user-supplied measurements with provenance, never distributed constants.
    for addition in observed_additions:
        if not addition.get("source"):
            raise ValueError("Supplemental peak requires an observation source")
        station = addition["station_id"]
        table = result[station]
        table = table.loc[~table.year.eq(addition["year"])].copy()
        result[station] = pd.concat(
            [
                table,
                pd.DataFrame(
                    [
                        dict(
                            year=addition["year"],
                            maximum=addition["water_m"],
                            coverage=np.nan,
                            coverage_JAS=np.nan,
                            used=True,
                            storm_quality_flag="Supplemental observed peak: " + addition["source"],
                        )
                    ]
                ),
            ],
            ignore_index=True,
        ).sort_values("year")
    return result, days, responses
