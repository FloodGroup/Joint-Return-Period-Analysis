"""Fit annual frequency curves and transfer a short tide record."""

from .annual_fit import Margin, annual_margin


def fit_annual(maxima, variable, settings):
    sample = maxima.loc[maxima.used & maxima.maximum.notna(), "maximum"].to_numpy(float)
    if variable == "R":
        sample = sample[sample > 0]
    if len(sample) < settings.minimum_annual_trial:
        return (
            dict(
                n_years=len(sample), diagnostic_only=True, reason="Fewer than eight annual maxima"
            ),
            [],
            [],
        )
    spec, candidates = annual_margin(
        sample, variable, seed=settings.annual_seed, reps=settings.annual_reps
    )
    if spec is None:
        return (
            dict(n_years=len(sample), diagnostic_only=True, reason="No eligible annual margin"),
            [],
            candidates,
        )
    diagnostic = bool(
        settings.synthetic or len(sample) < settings.minimum_annual_accepted or not spec["KS_pass"]
    )
    meta = dict(
        n_years=len(sample),
        spec=spec,
        diagnostic_only=diagnostic,
        used_years=maxima.loc[maxima.used, "year"].astype(int).tolist(),
    )
    rows = [
        dict(
            return_period_years=t,
            level=float(Margin(spec).ppf(1 - 1 / t)),
            diagnostic_only=diagnostic,
        )
        for t in settings.periods
    ]
    return meta, rows, candidates


def transfer_curve(short_maxima, long_maxima, long_design, minimum_year=None):
    short = short_maxima.loc[short_maxima.used].set_index("year").maximum
    long = long_maxima.loc[long_maxima.used].set_index("year").maximum
    common = short.index.intersection(long.index)
    if minimum_year is not None:
        common = common[common >= minimum_year]
    if len(common) < 3:
        raise ValueError("At least three concurrent annual maxima are required for transfer")
    delta = float((short.loc[common] - long.loc[common]).median())
    rows = [
        {**r, "level": r["level"] + delta, "transfer_delta_m": delta, "common_years": len(common)}
        for r in long_design
    ]
    return rows, dict(delta_m=delta, n_common=len(common), common_years=common.astype(int).tolist())
