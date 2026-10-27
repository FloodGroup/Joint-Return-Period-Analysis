"""Annual-maximum frequency analysis with numerical eligibility checks."""

import warnings
import numpy as np
from scipy import stats
from scipy.special import gammaln

POOL = [
    "gamma",
    "lognorm",
    "weibull_min",
    "expon",
    "genextreme",
    "genpareto",
    "logpearson3",
    "fisk",
    "invgauss",
    "gengamma",
    "pearson3",
    "burr12",
]
POSITIVE = {
    "gamma",
    "lognorm",
    "weibull_min",
    "expon",
    "genpareto",
    "fisk",
    "invgauss",
    "gengamma",
    "burr12",
}
ALPHA = 0.05


class Margin:
    def __init__(self, spec):
        self.spec = spec
        self.name = spec["distribution"]
        self.offset = float(spec.get("offset", 0))
        dn = "pearson3" if self.name == "logpearson3" else self.name
        self.dist = getattr(stats, dn)(*spec["parameters"])

    def base_cdf(self, x):
        x = np.asarray(x, float)
        if self.name == "logpearson3":
            with np.errstate(all="ignore"):
                return np.where(x > 0, self.dist.cdf(np.log(np.where(x > 0, x, 1))), 0.0)
        return self.dist.cdf(x)

    def cdf(self, x):
        return self.base_cdf(np.asarray(x, float) - self.offset)

    def ppf(self, u):
        q = np.clip(np.asarray(u, float), 1e-12, 1 - 1e-12)
        x = np.exp(self.dist.ppf(q)) if self.name == "logpearson3" else self.dist.ppf(q)
        return x + self.offset

    def logpdf(self, x):
        y = np.asarray(x, float) - self.offset
        with np.errstate(all="ignore"):
            if self.name == "logpearson3":
                return np.where(
                    y > 0,
                    self.dist.logpdf(np.log(np.where(y > 0, y, 1))) - np.log(np.where(y > 0, y, 1)),
                    -np.inf,
                )
            return self.dist.logpdf(y)


def validate_mlr(s):
    name = s["distribution"]
    p = s.get("parameters", [])
    if name in ["pearson3", "logpearson3"] and abs(p[0]) > 2:
        raise ValueError("Pearson III shape ineligible")
    if name in ["gamma", "weibull_min"] and p[0] < 1:
        raise ValueError("shape < 1")
    if name == "gengamma" and p[1] > 0 and (p[0] * p[1] < 1):
        raise ValueError("gengamma endpoint")
    if name in ["fisk", "burr12"] and p[0] < 1:
        raise ValueError("shape < 1")
    if name == "genpareto" and p[0] < -1:
        raise ValueError("GPD shape")
    if name == "genextreme" and p[0] > 1:
        raise ValueError("GEV shape")


def degenerate(s):
    """Reject nonfinite or weakly identified distribution parameters."""
    p = np.asarray(s.get("parameters", []), float)
    return bool(np.any(~np.isfinite(p)) or np.any(np.abs(p) > 10000.0))


def fit_spec(x, name, variable, offset=0.0, check=True):
    x = np.asarray(x, float)
    s = dict(distribution=name, variable=variable, offset=float(offset))
    y = x - offset
    d = getattr(stats, "pearson3" if name == "logpearson3" else name)
    if name == "logpearson3":
        if np.any(y <= 0):
            raise ValueError("nonpositive")
        yy = np.log(y)
        fixed = {}
    else:
        yy = y
        fixed = {"floc": 0.0} if variable == "R" and name in POSITIVE else {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        par = d.fit(yy, **fixed)
    s.update(
        parameters=list(map(float, par)),
        known_parameters={"loc": 0.0} if fixed else {},
        k=len(par) - len(fixed),
    )
    if check:
        validate_mlr(s)
        if degenerate(s):
            raise ValueError("degenerate parameters")
    m = Margin(s)
    lp = m.logpdf(x)
    if not np.isfinite(lp).all():
        raise ValueError("nonfinite likelihood")
    if check and variable == "R" and (m.base_cdf(0) > 1e-06):
        raise ValueError("negative rainfall probability")
    n = len(x)
    k = s["k"]
    s["loglik"] = float(lp.sum())
    den = n - k - 1
    s["AICc"] = 2 * k - 2 * s["loglik"] + (2 * k * (k + 1) / den if den > 0 else 0)
    return s


def gof(x, spec, reps, seed):
    name = spec["distribution"]
    sample = np.asarray(x) - spec.get("offset", 0)
    d = getattr(stats, "pearson3" if name == "logpearson3" else name)
    xx = np.log(sample) if name == "logpearson3" else sample
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = stats.goodness_of_fit(
            d,
            xx,
            known_params=spec.get("known_parameters", {}),
            statistic="ks",
            n_mc_samples=reps,
            random_state=np.random.default_rng(seed),
        )
    return dict(KS_D=float(res.statistic), KS_p=float(res.pvalue))


def sane_level(spec, variable):
    m = Margin(spec)
    try:
        q10 = float(m.ppf(1 - 1 / 10))
        q200 = float(m.ppf(1 - 1 / 200))
    except Exception:
        return False
    if not np.isfinite(q10) or not np.isfinite(q200):
        return False
    if variable == "H":
        return -1.0 <= q10 <= 5.5 and q10 <= q200 <= 8.0
    return 20 <= q10 <= 800 and q10 <= q200 <= 1500


def annual_margin(values, variable, seed=20261015, reps=99, exclude=()):
    """Select in AICc order with refitted KS and 10/200-year plausibility checks;
    otherwise first sane candidate (flagged 'diagnostic_tail')."""
    values = np.asarray(values, float)
    ranked, rows = ([], [])
    for i, name in enumerate([p for p in POOL if p not in exclude]):
        row = dict(variable=variable, distribution=name, selected=False)
        try:
            spec = fit_spec(values, name, variable, 0.0)
            row.update(AICc=spec["AICc"], status="fitted", parameters=spec["parameters"])
            row.update(gof(values, spec, reps, seed + i))
            spec.update(KS_D=row["KS_D"], KS_p=row["KS_p"], KS_pass=bool(row["KS_p"] >= ALPHA))
            ranked.append((spec, row))
        except Exception as exc:
            row.update(status="ineligible", reason=str(exc))
        rows.append(row)
    ranked.sort(key=lambda it: it[1].get("AICc", 1e99))
    chosen = None
    for spec, row in ranked:
        if row.get("KS_p", 0) >= ALPHA and sane_level(spec, variable):
            chosen = spec
            row["selected"] = True
            break
    if chosen is None:
        for spec, row in ranked:
            if sane_level(spec, variable):
                chosen = spec
                row["selected"] = True
                row["status"] = "diagnostic_tail"
                break
    return (chosen, rows)


def lmom_gev(x, T):
    """GEV by L-moments (Hosking 1990); returns quantiles for the return periods T (reference only)."""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    j = np.arange(1, n + 1)
    b0 = x.mean()
    b1 = ((j - 1) / (n - 1) * x).mean()
    b2 = ((j - 1) * (j - 2) / ((n - 1) * (n - 2)) * x).mean()
    l1, l2, l3 = (b0, 2 * b1 - b0, 6 * b2 - 6 * b1 + b0)
    t3 = l3 / l2
    c = 2 / (3 + t3) - np.log(2) / np.log(3)
    k = 7.859 * c + 2.9554 * c**2
    alpha = l2 * k / ((1 - 2 ** (-k)) * np.exp(gammaln(1 + k)))
    xi = l1 - alpha * (1 - np.exp(gammaln(1 + k))) / k
    T = np.atleast_1d(np.asarray(T, float))
    F = 1 - 1 / T
    return (xi + alpha / k * (1 - (-np.log(F)) ** k), dict(xi=xi, alpha=alpha, k=k))


def gumbel_mom(x, T):
    x = np.asarray(x, float)
    T = np.atleast_1d(np.asarray(T, float))
    a = np.sqrt(6) * x.std(ddof=1) / np.pi
    u = x.mean() - 0.5772 * a
    return u - a * np.log(-np.log(1 - 1 / T))
