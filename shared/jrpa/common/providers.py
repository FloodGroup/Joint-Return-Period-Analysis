"""Provider-export normalization with explicit, source-specific quality rules."""

from pathlib import Path
import numpy as np
import pandas as pd

POLICIES = {"jma_rain", "mlit_rain", "sf_rain", "jma_tide", "gesla_tide", "jodc_tide", "verified"}


def normalize_export(
    path,
    policy,
    variable,
    time_column="time_utc",
    value_column=None,
    flag_column="flag",
    scale=1.0,
    daily=False,
):
    """Screen tabular provider exports without interpolation or datum conversion."""
    if policy not in POLICIES or variable not in {"rain", "water"}:
        raise ValueError("Explicit provider quality policy and variable required")
    if (policy.endswith("rain") and variable != "rain") or (
        policy.endswith("tide") and variable != "water"
    ):
        raise ValueError("Provider quality policy does not match the variable")
    frame = pd.read_csv(path, low_memory=False)
    value_column = value_column or ("precip_mm" if variable == "rain" else "water_level_m")
    if not {time_column, value_column} <= set(frame):
        raise ValueError("Provider export is missing configured time or value columns")
    values = pd.to_numeric(frame[value_column], errors="coerce") * scale
    if policy == "verified":
        accepted = pd.Series(True, index=frame.index)
    else:
        if flag_column not in frame:
            raise ValueError("This provider policy requires an explicit quality-code column")
        flag = frame[flag_column].fillna("").astype(str).str.strip()
        accepted = {
            "jma_rain": lambda: flag.str[:1].isin(["8", "5"]),
            "mlit_rain": lambda: flag.isin(["", "*"]),
            "sf_rain": lambda: ~flag.str.contains("A", regex=False),
            "jma_tide": lambda: flag.isin(["", "R"]),
            "gesla_tide": lambda: flag.isin(["1", "1.0"]),
            "jodc_tide": lambda: ~flag.str.contains("R", regex=False),
        }[policy]()
    accepted &= np.isfinite(values)
    if variable == "rain":
        accepted &= values.ge(0)
    time = pd.to_datetime(frame[time_column], utc=not daily, format="mixed", errors="coerce")
    valid_time = time.notna()
    if not daily:
        valid_time &= time.dt.minute.eq(0) & time.dt.second.eq(0)
    values = values.where(accepted)
    value_key = "rain_mm" if variable == "rain" else "water_m"
    time_key = "date" if daily else "time_utc"
    result = pd.DataFrame({time_key: time, value_key: values}).loc[valid_time].copy()
    result["finite"] = result[value_key].notna()
    duplicates = int(result[time_key].duplicated().sum())
    result = result.sort_values(
        [time_key, "finite"], ascending=[True, False], kind="stable"
    ).drop_duplicates(time_key)
    audit = dict(
        provider_policy=policy,
        input_file=Path(path).name,
        input_rows=len(frame),
        rejected_quality_or_value=int((~accepted).sum()),
        rejected_time=int((~valid_time).sum()),
        duplicate_rows=duplicates,
        accepted_observations=int(result[value_key].notna().sum()),
        scale_factor=scale,
        datum_conversion_applied=False,
    )
    return result.drop(columns="finite").reset_index(drop=True), audit
