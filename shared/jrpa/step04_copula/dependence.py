"""IFM dependence estimation with refitted joint-CDF bootstrap tests."""

import json
import numpy as np
import pyvinecopulib as pv
from . import daily, marginals

FAMILIES = (
    "indep",
    "gaussian",
    "student",
    "clayton",
    "gumbel",
    "frank",
    "joe",
    "bb1",
    "bb6",
    "bb7",
    "bb8",
    "tawn",
)


def family_options():
    for family in FAMILIES:
        for rotation in (
            [0] if family in {"indep", "gaussian", "student", "frank"} else [0, 90, 180, 270]
        ):
            yield family, rotation


def make_u(x, specs):
    return np.asfortranarray(
        np.clip(
            np.column_stack(
                [marginals.Margin(specs[v]).cdf(x[:, j]) for j, v in enumerate(("R", "H"))]
            ),
            1e-10,
            1 - 1e-10,
        )
    )


def cvm(u, model):
    empirical = np.mean(np.all(u[None, :, :] <= u[:, None, :] + 1e-13, axis=2), axis=1)
    return float(len(u) * np.mean((empirical - model.cdf(u)) ** 2))


def select_margin(x, variable, offset, settings, seed):
    if settings.event_clock == "daily" and not settings.synthetic:
        return daily.select_margin(x, variable, offset, seed)
    pool = (
        marginals.POOL
        if variable == "R"
        else [n for n in marginals.POOL if n != "lognormal_mixture2"]
    )
    rows = []
    fitted = {}
    for i, name in enumerate(pool):
        row = dict(variable=variable, distribution=name, selected=False)
        try:
            spec = marginals.fit_spec(x, name, variable, offset, seed=seed + i)
            row.update({k: spec[k] for k in ("loglik", "k", "AIC", "AICc", "BIC")})
            row.update(marginals.gof(x, spec, settings.marginal_screen_reps, seed + i))
            row.update(status="fitted")
            fitted[name] = spec
        except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
            row.update(status="ineligible_or_failed", reason=str(exc))
        rows.append(row)
    ranked = sorted([r for r in rows if r["status"] == "fitted"], key=lambda r: r["AICc"])
    if not ranked:
        raise ValueError(f"No eligible {variable} margin")
    chosen = None
    for row in ranked:
        if row["KS_p"] < 0.02:
            continue
        try:
            row.update(
                marginals.gof(
                    x,
                    fitted[row["distribution"]],
                    settings.marginal_final_reps,
                    seed + 1000 + pool.index(row["distribution"]),
                )
            )
        except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
            row.update(KS_p=0.0, reason=str(exc))
        if row["KS_p"] >= 0.05:
            chosen = row
            break
    if chosen is None:
        chosen = ranked[0]
        if chosen.get("KS_replicates") != settings.marginal_final_reps:
            chosen.update(
                marginals.gof(
                    x, fitted[chosen["distribution"]], settings.marginal_final_reps, seed + 2000
                )
            )
    chosen["selected"] = True
    spec = fitted[chosen["distribution"]]
    spec.update(KS_D=chosen["KS_D"], KS_p=chosen["KS_p"], KS_pass=chosen["KS_p"] >= 0.05)
    return spec, rows


def joint_test(u, model, specs, settings, seed):
    statistic = cvm(u, model)
    trials = []
    kernel = daily if settings.event_clock == "daily" and not settings.synthetic else marginals
    for b in range(settings.joint_reps):
        try:
            sim = model.simulate(len(u), seeds=[seed + b])
            x = np.column_stack(
                [marginals.Margin(specs[v]).ppf(sim[:, j]) for j, v in enumerate(("R", "H"))]
            )
            estimated = {}
            for j, v in enumerate(("R", "H")):
                kwargs = {"check": False} if kernel is daily else {"check_eligibility": False}
                estimated[v] = kernel.fit_spec(
                    x[:, j],
                    specs[v]["distribution"],
                    v,
                    specs[v].get("offset", 0),
                    seed=seed + b,
                    **kwargs,
                )
            uu = make_u(x, estimated)
            fitted = pv.Bicop(family=model.family, rotation=model.rotation)
            fitted.fit(uu, controls=pv.FitControlsBicop(parametric_method="mle", num_threads=1))
            trials.append(dict(replicate=b, status="valid", CvM=cvm(uu, fitted)))
        except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
            trials.append(dict(replicate=b, status="failed", reason=str(exc)))
    valid = [r["CvM"] for r in trials if r["status"] == "valid"]
    enough = len(valid) >= 0.9 * settings.joint_reps
    p = (1 + sum(x >= statistic for x in valid)) / (1 + len(valid))
    return dict(
        CvM_statistic=statistic,
        CvM_p=p,
        bootstrap_valid=len(valid),
        bootstrap_replicates=settings.joint_reps,
        GOF_pass=bool(enough and p >= 0.05),
    ), trials


def fit_dependence(events, meta, settings, unit_number=0, water_override=False):
    x = events[["rain_mm", "water_max_m"]].to_numpy(float)
    if len(x) < 8 or not np.isfinite(x).all():
        raise ValueError("At least eight finite event pairs are required for a trial fit")
    is_daily = settings.event_clock == "daily"
    specs = {}
    margin_rows = []
    for j, v in enumerate(("R", "H")):
        seed = settings.seed + j * 100 + (0 if is_daily else unit_number * 10000)
        specs[v], rows = select_margin(
            x[:, j], v, meta["R_threshold_mm"] if v == "R" else 0.0, settings, seed
        )
        margin_rows.extend(rows)
    if water_override:
        from . import annual_fit

        ranked = []
        for i, name in enumerate(annual_fit.POOL):
            try:
                spec = annual_fit.fit_spec(x[:, 1], name, "H", 0.0)
                spec.update(annual_fit.gof(x[:, 1], spec, 49, settings.seed + 100 + i))
                ranked.append(spec)
            except (ValueError, RuntimeError, FloatingPointError, OverflowError):
                continue
        ranked.sort(key=lambda s: s["AICc"])
        if not ranked:
            raise ValueError("No nondegenerate water-level margin")
        selected = next((s for s in ranked if s["KS_p"] >= 0.05), ranked[0])
        selected["KS_pass"] = selected["KS_p"] >= 0.05
        specs["H"] = selected
        margin_rows = [r for r in margin_rows if r["variable"] != "H"] + [
            dict(**s, selected=s is selected) for s in ranked
        ]
    u = make_u(x, specs)
    candidates = []
    failures = []
    options = [("gaussian", 0)] if water_override else family_options()
    for family, rotation in options:
        row = dict(family=family, rotation=rotation, selected=False)
        try:
            model = pv.Bicop(family=getattr(pv.BicopFamily, family), rotation=rotation)
            model.fit(u, controls=pv.FitControlsBicop(parametric_method="mle", num_threads=1))
            ll = float(model.loglik(u))
            k = float(model.npars)
            n = len(u)
            row.update(
                AIC=2 * k - 2 * ll,
                AICc=2 * k - 2 * ll + 2 * k * (k + 1) / max(n - k - 1, 1),
                BIC=k * np.log(n) - 2 * ll,
                model_json=model.to_json(),
                status="fitted",
            )
            candidates.append(row)
        except (ValueError, RuntimeError) as exc:
            failures.append({**row, "status": "failed", "reason": str(exc)})
    candidates.sort(key=lambda r: r["AICc"])
    if not candidates:
        raise ValueError("No fitted copula candidates")
    tests = []
    all_trials = []
    winner = None
    for rank, row in enumerate(candidates, 1):
        model = pv.Bicop.from_json(row["model_json"])
        seed = settings.seed + (rank * 1000 if is_daily else unit_number * 1000)
        result, trials = joint_test(u, model, specs, settings, seed)
        tests.append(dict(rank=rank, family=row["family"], rotation=row["rotation"], **result))
        all_trials.extend(
            {**t, "rank": rank, "family": row["family"], "rotation": row["rotation"]}
            for t in trials
        )
        if result["GOF_pass"]:
            winner = row
            break
    winner = winner or candidates[0]
    winner["selected"] = True
    test = next(
        t for t in tests if (t["family"], t["rotation"]) == (winner["family"], winner["rotation"])
    )
    margin_pass = all(s["KS_pass"] for s in specs.values())
    record_pass = (
        len(x) >= settings.minimum_events
        and meta["effective_years"] >= settings.minimum_effective_years
    )
    fitted = {
        **meta,
        "margins": specs,
        "copula": json.loads(winner["model_json"]),
        "copula_family": winner["family"],
        "rotation": winner["rotation"],
        "AICc": winner["AICc"],
        "event_count": len(x),
        "marginal_KS_pass": margin_pass,
        "record_length_gate_pass": record_pass,
        **test,
        "diagnostic_only": bool(
            settings.synthetic or not (margin_pass and record_pass and test["GOF_pass"])
        ),
        "synthetic": settings.synthetic,
        "datum": settings.datum,
        "settings": settings.as_dict(),
    }
    return fitted, dict(
        marginal_candidates=margin_rows,
        copula_candidates=candidates + failures,
        copula_tests=tests,
        bootstrap_trials=all_trials,
    )
