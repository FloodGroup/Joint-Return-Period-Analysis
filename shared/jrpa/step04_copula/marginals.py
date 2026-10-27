"""Maximum-likelihood margins and refitted parametric-bootstrap tests."""

import warnings
import numpy as np
from scipy import stats
from scipy.special import logsumexp
from sklearn.mixture import GaussianMixture

SEED = 202610052
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


def fit_spec(x, name, variable, offset=0, mode="full", seed=SEED, check_eligibility=True):
    x = np.asarray(x, float)
    s = dict(distribution=name, variable=variable, offset=float(offset), mode="full", zero_mass=0.0)
    y = x - offset
    if not np.isfinite(y).all():
        raise ValueError("Nonfinite observations")
    if name == "lognormal_mixture2":
        if len(y) < 35 or np.any(y <= 0):
            raise ValueError("Mixture requires >=35 positive observations")
        g = GaussianMixture(
            2, covariance_type="full", n_init=8, reg_covar=1e-05, random_state=seed, max_iter=500
        )
        g.fit(np.log(y)[:, None])
        order = np.argsort(g.means_.ravel())
        s.update(
            weights=g.weights_[order].tolist(),
            log_means=g.means_.ravel()[order].tolist(),
            log_sd=np.sqrt(g.covariances_.ravel()[order]).tolist(),
            library="sklearn.mixture.GaussianMixture",
            k=5,
        )
        if not g.converged_ or (
            check_eligibility and (min(g.weights_) * len(y) < 8 or min(g.weights_) < 0.08)
        ):
            raise ValueError(
                "Mixture has nonconvergence or fewer than 8 effective observations in a component"
            )
        if check_eligibility and min(s["log_sd"]) < 0.05:
            raise ValueError("Near-collapsed mixture component")
    else:
        d = getattr(stats, "pearson3" if name == "logpearson3" else name)
        if name == "logpearson3":
            if np.any(y <= 0):
                raise ValueError("Log-Pearson III requires positive observations")
            yy = np.log(y)
            fixed = {}
        else:
            yy = y
            fixed = {"floc": 0.0} if variable in ["R", "Q"] and name in POSITIVE else {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            par = d.fit(yy, **fixed)
        s.update(
            parameters=list(map(float, par)),
            known_parameters={"loc": 0.0} if fixed else {},
            k=len(par) - len(fixed),
            library="scipy.stats." + d.name,
        )
        fr = d(*par)
        lo, hi = fr.support()
        if name != "logpearson3" and (not fixed):
            tolerance = 1e-06 * max(np.ptp(yy), 1.0)
            for b, ob in [(lo, min(yy)), (hi, max(yy))]:
                if np.isfinite(b) and abs(b - ob) < tolerance:
                    with np.errstate(all="ignore"):
                        bound_density = fr.pdf(b)
                    if check_eligibility and (not np.isfinite(bound_density)):
                        raise ValueError("Singular fitted endpoint at observed extremum")
    m = Margin(s)
    lp = m.logpdf(x)
    if not np.isfinite(lp).all():
        raise ValueError("Nonfinite fitted likelihood")
    if check_eligibility and variable in ["R", "Q"] and (m.base_cdf(0) > 1e-06):
        raise ValueError("Fit assigns >1e-6 probability to negative rainfall/flow")
    s["loglik"] = float(lp.sum())
    n = len(x)
    k = s["k"]
    s.update(
        n=n,
        AIC=2 * k - 2 * s["loglik"],
        AICc=2 * k - 2 * s["loglik"] + 2 * k * (k + 1) / (n - k - 1),
        BIC=k * np.log(n) - 2 * s["loglik"],
    )
    s["mode"] = mode
    if check_eligibility:
        validate_mlr_margin(s, x)
    return s


def validate_mlr_margin(s, x):
    """Exclude mathematically unbounded endpoint modes and unidentified shifts.

    These are design-density eligibility checks, not KS rejections. SciPy's
    Pearson III generic support reports (-inf,inf); its actual gamma endpoint
    is obtained from the documented parameterization, without fitting a model.
    """
    name = s["distribution"]
    p = s.get("parameters", [])
    if name in ["pearson3", "logpearson3"] and abs(p[0]) > 2:
        raise ValueError(
            "MLR ineligible: Pearson III shape 4/skew^2 <1 gives an unbounded finite-endpoint density"
        )
    if name in ["gamma", "weibull_min"] and p[0] < 1:
        raise ValueError("MLR ineligible: unbounded lower-endpoint density (shape <1)")
    if name == "gengamma" and p[1] > 0 and (p[0] * p[1] < 1):
        raise ValueError("MLR ineligible: unbounded lower-endpoint density (a*c <1)")
    if name in ["fisk", "burr12"] and p[0] < 1:
        raise ValueError("MLR ineligible: unbounded lower-endpoint density (shape <1)")
    if name == "genpareto" and p[0] < -1:
        raise ValueError("MLR ineligible: unbounded upper-endpoint density (GPD shape <-1)")
    if name == "genextreme" and p[0] > 1:
        raise ValueError("MLR ineligible: unbounded upper-endpoint density (SciPy GEV shape >1)")
    if name in POSITIVE and (not s.get("known_parameters")) and (len(p) > 1):
        span = max(float(np.ptp(x)), 1e-08)
        if abs(p[-2] - float(np.min(x))) > 100 * span:
            raise ValueError(
                "MLR ineligible: free location >100 observed ranges from observations; weakly identified shift"
            )


def ad_stat(cdf):
    u = np.clip(np.sort(cdf), 1e-12, 1 - 1e-12)
    n = len(u)
    i = np.arange(1, n + 1)
    return float(-n - np.sum((2 * i - 1) * (np.log(u) + np.log1p(-u[::-1]))) / n)


def gof(x, spec, reps, seed):
    """Refitted KS using SciPy; mixtures use an explicit refitted library bootstrap."""
    m = Margin(spec)
    name = spec["distribution"]
    sample = np.asarray(x) - spec.get("offset", 0)
    observedcdf = m.base_cdf(sample)
    D = float(stats.kstest(observedcdf, "uniform").statistic)
    if name != "lognormal_mixture2":
        d = getattr(stats, "pearson3" if name == "logpearson3" else name)
        xx = np.log(sample) if name == "logpearson3" else sample
        names = (d.shapes.split(", ") if d.shapes else []) + ["loc", "scale"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = stats.goodness_of_fit(
                d,
                xx,
                known_params=spec.get("known_parameters", {}),
                fit_params=dict(zip(names, spec["parameters"])),
                statistic="ks",
                n_mc_samples=reps,
                rng=np.random.default_rng(seed),
            )
        return dict(
            KS_D=float(res.statistic),
            KS_p=float(res.pvalue),
            KS_replicates=reps,
            KS_critical95=float(np.quantile(res.null_distribution, 0.95)),
            AD_statistic=ad_stat(observedcdf),
            KS_method="SciPy goodness_of_fit; parameters refitted in every bootstrap",
        )
    rng = np.random.default_rng(seed)
    null = []
    adnull = []
    fail = 0
    for b in range(reps):
        try:
            xb = m.ppf(rng.uniform(size=len(x)))
            sb = fit_spec(
                xb,
                name,
                spec["variable"],
                spec.get("offset", 0),
                seed=seed + b + 1,
                check_eligibility=False,
            )
            ub = Margin(sb).base_cdf(xb - spec.get("offset", 0))
            null.append(float(stats.kstest(ub, "uniform").statistic))
            adnull.append(ad_stat(ub))
        except Exception:
            fail += 1
    if len(null) < 0.95 * reps:
        raise ValueError("Mixture bootstrap has >5% failures")
    return dict(
        KS_D=D,
        KS_p=(1 + sum((t >= D for t in null))) / (1 + len(null)),
        KS_replicates=len(null),
        KS_critical95=float(np.quantile(null, 0.95)),
        AD_statistic=ad_stat(observedcdf),
        AD_p=(1 + sum((t >= ad_stat(observedcdf) for t in adnull))) / (1 + len(adnull)),
        AD_replicates=len(adnull),
        KS_method="Explicit refitted bootstrap; every mixture refitted by sklearn",
        bootstrap_failures=fail,
    )
