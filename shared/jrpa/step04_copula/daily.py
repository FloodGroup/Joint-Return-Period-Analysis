"""Daily-event marginal estimation for the Greater Bay Area."""

import warnings
import numpy as np
from scipy import stats
from scipy.special import logsumexp
from sklearn.mixture import GaussianMixture

SEED = 20261011
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
    "lognormal_mixture2",
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
        self.p0 = float(spec.get("zero_mass", 0))
        if self.name == "lognormal_mixture2":
            self.weights = np.asarray(spec["weights"])
            self.means = np.asarray(spec["log_means"])
            self.sd = np.asarray(spec["log_sd"])
            logx = np.linspace(
                np.min(self.means - 12 * self.sd), np.max(self.means + 12 * self.sd), 20000
            )
            f = np.sum(
                self.weights * stats.norm.cdf((logx[:, None] - self.means) / self.sd), axis=1
            )
            unique = np.r_[True, np.diff(f) > 0]
            self.inverse = (f[unique], logx[unique])
        else:
            dist_name = "pearson3" if self.name == "logpearson3" else self.name
            self.dist = getattr(stats, dist_name)(*spec["parameters"])

    def base_cdf(self, x):
        x = np.asarray(x, float)
        if self.name == "lognormal_mixture2":
            with np.errstate(all="ignore"):
                return np.where(
                    x > 0,
                    np.sum(
                        self.weights
                        * stats.norm.cdf((np.log(x)[..., None] - self.means) / self.sd),
                        axis=-1,
                    ),
                    0.0,
                )
        if self.name == "logpearson3":
            with np.errstate(all="ignore"):
                return np.where(x > 0, self.dist.cdf(np.log(x)), 0.0)
        return self.dist.cdf(x)

    def cdf(self, x):
        x = np.asarray(x, float)
        return self.p0 + (1 - self.p0) * self.base_cdf(x - self.offset)

    def ppf(self, u):
        u = np.asarray(u, float)
        q = np.clip((u - self.p0) / (1 - self.p0), 1e-12, 1 - 1e-12)
        if self.name == "lognormal_mixture2":
            x = np.exp(np.interp(q, *self.inverse))
        elif self.name == "logpearson3":
            x = np.exp(self.dist.ppf(q))
        else:
            x = self.dist.ppf(q)
        return x + self.offset

    def logpdf(self, x):
        x = np.asarray(x, float)
        y = x - self.offset
        with np.errstate(all="ignore"):
            if self.name == "lognormal_mixture2":
                lp = logsumexp(
                    np.log(self.weights)
                    + stats.norm.logpdf(np.log(y)[..., None], self.means, self.sd),
                    axis=-1,
                ) - np.log(y)
                lp = np.where(y > 0, lp, -np.inf)
            elif self.name == "logpearson3":
                lp = np.where(y > 0, self.dist.logpdf(np.log(y)) - np.log(y), -np.inf)
            else:
                lp = self.dist.logpdf(y)
        return lp + np.log(1 - self.p0)


def validate_mlr(s):
    name = s["distribution"]
    p = s.get("parameters", [])
    if name in ["pearson3", "logpearson3"] and abs(p[0]) > 2:
        raise ValueError("Pearson III shape ineligible")
    if name in ["gamma", "weibull_min"] and p[0] < 1:
        raise ValueError("shape < 1")
    if name == "gengamma" and p[1] > 0 and p[0] * p[1] < 1:
        raise ValueError("gengamma endpoint")
    if name in ["fisk", "burr12"] and p[0] < 1:
        raise ValueError("shape < 1")
    if name == "genpareto" and p[0] < -1:
        raise ValueError("GPD shape")
    if name == "genextreme" and p[0] > 1:
        raise ValueError("GEV shape")


def fit_spec(x, name, variable, offset=0, seed=SEED, check=True):
    x = np.asarray(x, float)
    s = dict(distribution=name, variable=variable, offset=float(offset), zero_mass=0.0)
    y = x - offset
    if name == "lognormal_mixture2":
        if len(y) < 35 or np.any(y <= 0):
            raise ValueError("mixture size")
        g = GaussianMixture(
            2, covariance_type="full", n_init=8, reg_covar=1e-5, random_state=seed, max_iter=500
        )
        g.fit(np.log(y)[:, None])
        order = np.argsort(g.means_.ravel())
        s.update(
            weights=g.weights_[order].tolist(),
            log_means=g.means_.ravel()[order].tolist(),
            log_sd=np.sqrt(g.covariances_.ravel()[order]).tolist(),
            k=5,
        )
        if not g.converged_ or (
            check
            and (min(g.weights_) * len(y) < 8 or min(g.weights_) < 0.08 or min(s["log_sd"]) < 0.05)
        ):
            raise ValueError("mixture collapsed")
    else:
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
    if check and name != "lognormal_mixture2":
        validate_mlr(s)
    m = Margin(s)
    lp = m.logpdf(x)
    if not np.isfinite(lp).all():
        raise ValueError("nonfinite likelihood")
    if check and variable == "R" and m.base_cdf(0) > 1e-6:
        raise ValueError("negative rainfall probability")
    n = len(x)
    k = s["k"]
    s["loglik"] = float(lp.sum())
    den = n - k - 1
    s["AICc"] = 2 * k - 2 * s["loglik"] + (2 * k * (k + 1) / den if den > 0 else 0)
    return s


def gof(x, spec, reps, seed):
    m = Margin(spec)
    name = spec["distribution"]
    sample = np.asarray(x) - spec.get("offset", 0)
    D = float(stats.kstest(m.base_cdf(sample), "uniform").statistic)
    if name == "lognormal_mixture2":
        return dict(KS_D=D, KS_p=np.nan, KS_replicates=0)
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
    return dict(KS_D=float(res.statistic), KS_p=float(res.pvalue), KS_replicates=reps)


def select_margin(x, variable, offset, seed):
    pool = POOL if variable == "R" else [n for n in POOL if n != "lognormal_mixture2"]
    rows = []
    fitted = {}
    for i, name in enumerate(pool):
        row = dict(variable=variable, distribution=name, selected=False)
        try:
            s = fit_spec(x, name, variable, offset, seed + i)
            row.update(AICc=s["AICc"], status="fitted")
            row.update(gof(x, s, 49, seed + i))
            fitted[name] = s
        except Exception as exc:
            row.update(status="ineligible", reason=str(exc))
        rows.append(row)
    ranked = sorted([r for r in rows if r["status"] == "fitted"], key=lambda r: r["AICc"])
    if not ranked:
        for i, name in enumerate(pool):
            try:
                s = fit_spec(x, name, variable, offset, seed + i, check=False)
                row = dict(
                    variable=variable,
                    distribution=name,
                    selected=False,
                    AICc=s["AICc"],
                    status="fitted_unchecked",
                )
                fitted[name] = s
                rows.append(row)
            except Exception:
                continue
        ranked = sorted([r for r in rows if "AICc" in r], key=lambda r: r["AICc"])
    if not ranked:
        reasons = [r.get("reason") for r in rows if r.get("status") == "ineligible"][:6]
        raise RuntimeError("no margin " + variable + " " + str(reasons))
    chosen = None
    for row in ranked:
        if row.get("KS_p", 0) >= ALPHA:
            chosen = row
            break
    if chosen is None:
        chosen = ranked[0]
    if "KS_p" not in chosen:
        chosen.update(gof(x, fitted[chosen["distribution"]], 49, seed + 17))
    chosen["selected"] = True
    s = fitted[chosen["distribution"]]
    s.update(KS_D=chosen["KS_D"], KS_p=chosen["KS_p"], KS_pass=bool(chosen.get("KS_p", 0) >= ALPHA))
    return s, rows
