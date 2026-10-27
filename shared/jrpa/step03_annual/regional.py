"""Distance-limited index-rainfall frequency analysis using regional L-moments."""

import numpy as np
from scipy.special import gammaln
from jrpa.gba.settings import (
    REGIONAL_RADIUS_KM,
    REGIONAL_H1_LIMIT,
    REGIONAL_REPLICATES,
    REGIONAL_SEED,
)


def distance_km(a, b):
    la1, lo1, la2, lo2 = map(np.radians, (*a, *b))
    h = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return float(6371 * 2 * np.arcsin(np.sqrt(h)))


def l_moments(values):
    x = np.sort(np.asarray(values, float))
    n = len(x)
    if n < 3 or not np.isfinite(x).all() or x.mean() <= 0:
        raise ValueError("L-moments require three finite observations with positive mean")
    j = np.arange(1, n + 1)
    b0 = x.mean()
    b1 = (((j - 1) / (n - 1)) * x).mean()
    b2 = (((j - 1) * (j - 2) / ((n - 1) * (n - 2))) * x).mean()
    l2 = 2 * b1 - b0
    if l2 <= 0:
        raise ValueError("Positive L-scale required")
    return float(b0), float(l2 / b0), float((6 * b2 - 6 * b1 + b0) / l2)


def growth_curve(tau, t3, periods):
    c = 2 / (3 + t3) - np.log(2) / np.log(3)
    k = 7.8590 * c + 2.9554 * c * c
    if abs(k) < 1e-7:
        a = tau / np.log(2)
        xi = 1 - a * np.euler_gamma
        q = {t: float(xi - a * np.log(-np.log(1 - 1 / t))) for t in periods}
    else:
        a = tau * k / ((1 - 2 ** (-k)) * np.exp(gammaln(1 + k)))
        xi = 1 - a * (1 - np.exp(gammaln(1 + k))) / k
        q = {t: float(xi + a / k * (1 - (-np.log(1 - 1 / t)) ** k)) for t in periods}
    return q, (xi, a, k)


def heterogeneity(donors, replicates=REGIONAL_REPLICATES, seed=REGIONAL_SEED):
    n = np.array([len(d["maxima"]) for d in donors])
    ratios = np.array([l_moments(d["maxima"])[1:] for d in donors])
    weights = n / n.sum()
    tau, t3 = weights @ ratios
    observed = np.sqrt(np.sum(n * (ratios[:, 0] - tau) ** 2) / n.sum())
    _, (xi, a, k) = growth_curve(tau, t3, (100,))
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(replicates):
        tc = []
        for size in n:
            u = rng.random(size)
            x = (
                (xi - a * np.log(-np.log(u)))
                if abs(k) < 1e-7
                else xi + a / k * (1 - (-np.log(u)) ** k)
            )
            tc.append(l_moments(x)[1])
        tc = np.asarray(tc)
        center = weights @ tc
        null.append(np.sqrt(np.sum(n * (tc - center) ** 2) / n.sum()))
    return float((observed - np.mean(null)) / np.std(null))


def regional_design(local_maxima, location, donors, periods=(10, 50, 100, 200), factor=1.13):
    candidates = []
    for donor in donors:
        distance = distance_km(location, donor["location"])
        if distance <= REGIONAL_RADIUS_KM:
            candidates.append({**donor, "distance_km": distance})
    candidates.sort(key=lambda d: d["distance_km"])
    if not candidates:
        raise ValueError("No eligible regional rainfall stations within 60 km")
    dropped = []
    h = heterogeneity(candidates) if len(candidates) > 1 else None
    while len(candidates) > 1 and h >= REGIONAL_H1_LIMIT:
        dropped.append(candidates.pop()["station_id"])
        h = heterogeneity(candidates) if len(candidates) > 1 else None
    lengths = np.array([len(d["maxima"]) for d in candidates])
    ratios = np.array([l_moments(d["maxima"])[1:] for d in candidates])
    tau, t3 = (lengths / lengths.sum()) @ ratios
    growth, _ = growth_curve(tau, t3, periods)
    local = np.asarray(local_maxima, float)
    if not len(local) or not np.isfinite(local).all():
        raise ValueError("Finite local annual maxima required")
    index = float(local.mean())
    values = [dict(return_period_years=t, level=index * growth[t] * factor) for t in periods]
    audit = dict(
        method="Regional GEV L-moments / index rainfall",
        index_mm=index,
        H1=h,
        L_CV=tau,
        L_skew=t3,
        station_years=int(lengths.sum()),
        members=[
            dict(station_id=d["station_id"], distance_km=d["distance_km"], n_years=len(d["maxima"]))
            for d in candidates
        ],
        dropped=dropped,
        growth=growth,
    )
    return values, audit
