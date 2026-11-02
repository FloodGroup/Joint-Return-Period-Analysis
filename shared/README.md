# Analysis workflow

Run the study-area entry script: [GBA](../GBA/run_analysis.py), [SanF](../SanF/run_analysis.py) or [Tokyo](../Tokyo/run_analysis.py). Each calls the five steps below. Regional settings and station assignments are stored alongside the entry script. [Methods](methods.md) gives the equations and detailed regional rules; [Inputs and sources](inputs.md) describes the private observation files.

| Step | Processing | Code |
|--:|:--|:--|
| 01 | Prepare observations and vertical references | [step01_prepare/run.py](jrpa/step01_prepare/run.py) |
| 02 | Merge rainfall and construct event samples | [step02_events/run.py](jrpa/step02_events/run.py) |
| 03 | Fit annual rainfall and water-level frequency curves | [step03_annual/run.py](jrpa/step03_annual/run.py) |
| 04 | Fit event margins and copula dependence | [step04_copula/run.py](jrpa/step04_copula/run.py) |
| 05 | Calculate design events and export results | [step05_export/run.py](jrpa/step05_export/run.py) |

## 01. Prepare observations

Read the observation manifest, validate timestamps and units, apply documented quality exclusions, and convert water levels to CN1985 (GBA), NAVD88 (SanF) or Tokyo Peil (Tokyo). Datum offsets have explicit validity intervals and evidence sources. Missing observations remain missing. This step also prepares tide annual samples, loads GBA regional rainfall records and records input hashes.

## 02. Construct event samples

Merge rainfall observations in the declared station-priority order and record each station's actual contribution. Prepare rainfall annual maxima separately from the paired event sample.

| Region | Rainfall–water-level pairing | Event selection |
|:--|:--|:--|
| SanF and Tokyo | Separate rainfall events by at least 12 dry hours. Use the complete total for events shorter than 24 hours; otherwise use the maximum complete 24-hour window. Pair with the highest water level in that window, requiring at least 75% water-level coverage. | Select strict exceedances of the 95th percentile of rainfall-complete event totals, calculated before water-coverage rejection. Cluster events with adjacent water peaks less than 72 hours apart and retain the largest rainfall. |
| GBA | Pair daily rainfall with the same day's highest water level. | Select strict exceedances of the 95th percentile of positive daily rainfall on the declared analysis clock. Cluster adjacent exceedance dates less than three days apart and retain the largest rainfall, using water level to resolve ties. |

Effective paired years are calculated from simultaneous valid observations. The selected-event rate is the number of retained events divided by this exposure.

## 03. Perform univariate return-period analysis

Perform separate univariate return-period analyses for rainfall (S1) and coastal water level (S2), using annual maxima screened independently for each variable. S2 is fitted once per tide station and reused by units sharing that record. The target return periods are 10, 50, 100 and 200 years.

The annual candidate pool contains **Gamma, Lognormal, Weibull, Exponential, Generalized Extreme Value (GEV), Generalized Pareto (GPD), Log-Pearson III, Log-logistic, Inverse Gaussian, Generalized Gamma, Pearson III and Burr XII**. Parameters are estimated by maximum likelihood; eligible candidates are ranked by the corrected Akaike information criterion (AICc) [1].

A refitted Kolmogorov–Smirnov (KS) test uses 99 Monte Carlo replicates [2]. Select the first candidate in AICc order that has p ≥ 0.05 and passes the 10- and 200-year plausibility checks. A candidate passing only the plausibility checks is retained for diagnosis.

For GBA mainland rainfall units, a regional L-moment procedure replaces the local distribution search: the local mean annual maximum scales a regional GEV growth curve [3]. Donors are restricted to 60 km, with iterative heterogeneity screening, and a factor of 1.13 converts daily to estimated moving-24-hour rainfall. Specified short tide records use donor curves adjusted by the median difference of concurrent annual maxima.

## 04. Fit event margins and copulas

Fit rainfall excess above its selection threshold and water level on the paired event sample. The event margins use the same 12 distribution families as Step 03; rainfall additionally admits a **two-component lognormal mixture**. Mixtures are fitted by expectation–maximization. Numerical and density-shape eligibility checks precede AICc ranking.

The copula pool contains **Independence, Gaussian, Student-t, Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7, BB8 and Tawn**. Clayton, Gumbel, Joe, BB1, BB6, BB7, BB8 and Tawn also use 90°, 180° and 270° rotations, giving 36 combinations. Estimate margins first, then copula parameters by maximum likelihood on the transformed observations using pyvinecopulib [5]. GBA 10a uses the annual-family eligibility rules for its event water-level margin and a Gaussian copula.

| Check | GBA | SanF and Tokyo |
|:--|:--|:--|
| Event-margin KS | 49 refitted replicates | 49 screening replicates; candidates with screening p ≥ 0.02 proceed to 999-replicate testing in AICc order |
| Joint-CDF Cramér–von Mises (CvM) | 49 bootstrap replicates | 99 bootstrap replicates |
| Formal compound estimate | Both marginal p values ≥ 0.05; joint p ≥ 0.05 with at least 90% valid replicates; at least 30 event pairs and 10 effective paired years | Same criteria |

The joint test compares the empirical and fitted joint CDFs, using a parametric bootstrap with both margins and the fixed copula family refitted in each replicate. Copulas are tested in AICc order; the first passing candidate is selected. The CvM and bootstrap foundations are discussed in [4]; the counts and thresholds above are repository settings.

GBA's lognormal mixture has no KS acceptance p value and cannot satisfy the formal marginal criterion. If no copula passes, the best-ranked fitted candidate remains available for diagnosis. Exported results retain their diagnostic status and reason.

## 05. Calculate and export designs

Calibrate a Kendall return-period contour using the selected-event rate and simulated copula probabilities [6]. On each contour, select the rainfall–water-level combination with the highest **joint density in physical units**, including both marginal densities. SanF and Tokyo use 1,048,576 calibration samples and an equally sized independent validation sample; GBA uses 262,144 calibration samples.

Combine S1, S2 and Kendall designs with the supplied climate inputs to define rainfall-dominated, storm-surge-dominated and compound rainfall–storm-surge scenarios. The marine variable is total coastal water level.

## References

1. Hurvich, C. M., & Tsai, C.-L. (1989). Regression and time series model selection in small samples. *Biometrika*, **76**(2), 297–307. [doi:10.1093/biomet/76.2.297](https://doi.org/10.1093/biomet/76.2.297).
2. SciPy developers. *SciPy 1.15.3: scipy.stats.goodness_of_fit*. [Refitted Monte Carlo goodness-of-fit procedure](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.goodness_of_fit.html).
3. Hosking, J. R. M., & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press. [Publisher's description](https://www.cambridge.org/core/books/abs/regional-frequency-analysis/regional-frequency-analysis/382FDB3E981D331A9A5EB11A1D33CEEC).
4. Genest, C., Rémillard, B., & Beaudoin, D. (2009). Goodness-of-fit tests for copulas: A review and a power study. *Insurance: Mathematics and Economics*, **44**(2), 199–213. [doi:10.1016/j.insmatheco.2007.10.005](https://doi.org/10.1016/j.insmatheco.2007.10.005).
5. pyvinecopulib developers. *pyvinecopulib 0.7.5: FitControlsBicop*. [Maximum-likelihood fitting controls](https://vinecopulib.github.io/pyvinecopulib/_generate/pyvinecopulib.FitControlsBicop.__init__.html).
6. Salvadori, G., De Michele, C., & Durante, F. (2011). On the return period and design in a multivariate framework. *Hydrology and Earth System Sciences*, **15**, 3293–3305. [doi:10.5194/hess-15-3293-2011](https://doi.org/10.5194/hess-15-3293-2011).

[Main runner](jrpa/main.py) executes the same steps for `jrpa run` and `jrpa demo`. Shared settings, file handling and quality controls are in `jrpa/common/`. Study-area maps and observation-coverage figures are in each region's `figures/` folder.
