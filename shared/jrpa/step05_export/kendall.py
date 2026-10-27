"""Ordinary Kendall calibration and maximum physical joint-density design."""

import json
import numpy as np
import pyvinecopulib as pv
from scipy.optimize import brentq, minimize_scalar
from jrpa.step04_copula.marginals import Margin


def design_points(meta, settings, unit_number=0):
    model = pv.Bicop.from_json(json.dumps(meta["copula"]))
    margins = [Margin(meta["margins"][v]) for v in ("R", "H")]
    rate = meta["events_per_year"]
    if not np.isfinite(rate) or rate <= 0:
        raise ValueError("Positive event rate required")
    seed = settings.seed + (0 if settings.event_clock == "daily" else unit_number * 100)
    scores = model.cdf(model.simulate(settings.kendall_samples, qrng=True, seeds=[seed + 11]))
    validation = (
        model.cdf(model.simulate(settings.kendall_samples, qrng=True, seeds=[seed + 71]))
        if settings.independent_calibration
        else None
    )
    rows = []
    curves = []
    search = []
    calibration = []
    for period in settings.periods:
        probability = 1 / (rate * period)
        if not 0 < probability < 1:
            raise ValueError("Return period is incompatible with the event rate")
        critical = float(np.quantile(scores, 1 - probability))
        achieved = float(np.mean(validation > critical)) if validation is not None else None
        error = abs(achieved / probability - 1) if achieved is not None else None
        if error is not None and error >= 0.05:
            raise ValueError(f"Independent Kendall probability error exceeds 5% at T={period}")

        def v_at(u):
            if model.cdf(np.asfortranarray([[u, 1 - 1e-11]]))[0] < critical:
                return np.nan
            bounds = (
                (1e-8, 1 - 1e-8)
                if settings.kendall_search == "daily_grid"
                else (max(critical, 1e-10), 1 - 1e-11)
            )
            return brentq(
                lambda v: float(model.cdf(np.asfortranarray([[u, v]]))[0]) - critical,
                *bounds,
                xtol=1e-10 if settings.kendall_search == "daily_grid" else 1e-12,
            )

        def logdensity(u):
            v = v_at(float(u))
            if not np.isfinite(v):
                return -np.inf
            cube = np.asfortranarray([[float(u), v]])
            x = [float(margins[j].ppf(cube[0, j])) for j in range(2)]
            return float(
                np.log(max(model.pdf(cube)[0], 1e-300))
                + sum(margins[j].logpdf(x[j]) for j in range(2))
            )

        if settings.kendall_search == "daily_grid":
            grid = np.linspace(min(0.999, max(critical + 1e-4, 0.02)), 0.999, 80)
        else:
            low = critical + 1e-9
            upper = 1 - 1e-9
            grid = np.sort(
                np.unique(
                    np.r_[
                        np.linspace(low, upper, 321),
                        1 - np.geomspace(1e-9, max(1 - critical - 1e-9, 2e-9), 321),
                    ]
                )
            )
            grid = grid[(grid >= low) & (grid <= upper)]
        density = np.array([logdensity(u) for u in grid])
        if not np.isfinite(density).any():
            raise ValueError("No feasible Kendall layer")
        best = int(np.nanargmax(density))
        if settings.kendall_search == "daily_grid":
            opt = minimize_scalar(
                lambda u: -logdensity(u),
                bounds=(grid[max(0, best - 1)], grid[min(len(grid) - 1, best + 1)]),
                method="bounded",
            )
            choices = [(-float(opt.fun), float(opt.x))]
        else:
            choices = [(density[best], grid[best])]
            for j in range(1, len(grid) - 1):
                if (
                    np.isfinite(density[j])
                    and density[j] >= density[j - 1]
                    and density[j] >= density[j + 1]
                ):
                    opt = minimize_scalar(
                        lambda u: -logdensity(u),
                        bounds=(grid[j - 1], grid[j + 1]),
                        method="bounded",
                        options={"xatol": 1e-13},
                    )
                    choices.append((-float(opt.fun), float(opt.x)))
        lp, u = max(choices)
        v = v_at(u)
        r = float(margins[0].ppf(u))
        h = float(margins[1].ppf(v))
        residual = abs(float(model.cdf(np.asfortranarray([[u, v]]))[0]) - critical)
        if residual > 1e-8 or r < meta["R_threshold_mm"] - 1e-7:
            raise ValueError("Kendall layer or rainfall-threshold residual failed")
        rows.append(
            dict(
                unit=meta["unit"],
                return_period_years=period,
                scenario="Kendall",
                R_mm=r * settings.rainfall_factor_24h,
                H_m=h,
                R_quantile=u,
                H_quantile=v,
                critical_CDF_score=critical,
                layer_residual=residual,
                log_joint_density=lp,
                diagnostic_only=meta["diagnostic_only"],
                synthetic=settings.synthetic,
            )
        )
        for u0 in grid[np.isfinite(density)]:
            curves.append(
                dict(
                    return_period_years=period,
                    R_mm=float(margins[0].ppf(u0)) * settings.rainfall_factor_24h,
                    H_m=float(margins[1].ppf(v_at(u0))),
                )
            )
        search.extend(dict(return_period_years=period, u=u0, log_density=d) for d, u0 in choices)
        calibration.append(
            dict(
                return_period_years=period,
                critical_CDF_score=critical,
                target_probability=probability,
                validation_probability=achieved,
                relative_error=error,
            )
        )
    return (
        rows,
        curves,
        dict(
            n_calibration=settings.kendall_samples,
            n_validation=settings.kendall_samples if validation is not None else 0,
            independent_validation=validation is not None,
            levels=calibration,
            search=search,
        ),
    )
