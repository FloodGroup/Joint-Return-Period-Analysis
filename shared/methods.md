# Methods

Each region's `run_analysis.py` is the program to run. It calls the five folders under `shared/jrpa/` in the order below. `jrpa.main.run` calls the same five functions. Each region's `settings.py` defines its parameters and station rules; `stations.json` supplies its station catalog. Rainfall and coastal water level are denoted by $R$ and $H$.

## 01. Prepare observations

`step01_prepare` reads the private manifest, checks the time axis and vertical reference, applies documented exclusions and assembles water-level annual samples. It also loads rainfall stations used for GBA regional frequency analysis. Input paths and hashes accompany the output.

Rainfall follows a fixed station-priority sequence. Valid zeros remain observations; missing values remain missing. Water levels must use the regional analysis datum before fitting. Epoch-specific conversions require a documented offset and date interval.

## 02. Construct event samples

`step02_events` merges rainfall sources, builds paired rainfall and water-level samples, calculates paired exposure and prepares rainfall annual maxima.

For San Francisco and Tokyo, positive accumulation intervals separated by at least 12 hours define rainfall events. An event shorter than 24 hours uses its complete event total and the highest coastal water level over the same window. Longer events use their maximum complete 24-hour rainfall window and its coincident highest water level. At least 75% of water-level hours must be observed in the selected window.

The rainfall threshold is the 95th percentile of rainfall-complete event totals before water-coverage rejection and declustering. Strict exceedances are clustered by adjacent water-peak separations below 72 hours. The event with the greatest rainfall is retained, with the earlier window resolving ties. An unobserved dry interval is not treated as a confirmed zero-rainfall interval.

For GBA, the threshold is the 95th percentile of positive daily rainfall within the declared analysis clock. Same-day daily highest water levels are paired with rainfall exceedances. Adjacent exceedance dates separated by fewer than three days form clusters. The largest rainfall is retained; water level resolves rainfall ties. Explicitly invalid selected pairs are removed without replacing them with another storm in the cluster.

$$
Y_{\mathrm{eff}}=\frac{N_{\mathrm{paired}}}{365.2425\,s},
\qquad \lambda=\frac{N_{\mathrm{events}}}{Y_{\mathrm{eff}}},
$$

where $s=24$ for hourly observations and $s=1$ for daily observations. Exposure counts simultaneous valid observations, not the span between first and last records.

## 03. Fit annual frequency curves

`step03_annual` estimates full-record S1 rainfall and S2 water-level designs. It fits each tide station once, applies the specified short-record transfers and assigns the resulting S2 curve to all units using that station.

S1 uses full-record rainfall annual maxima, independently of water-level availability. San Francisco and Tokyo use complete sliding 24-hour totals; years follow fixed PST and JST October-to-September water years, respectively. S2 is fitted once per tide station on full-record annual maxima. Coverage is variable-specific. At least eight maxima are needed for a trial fit; fewer than 15 or a failed KS test implies diagnostic status.

Annual distribution selection uses the 12 single-family candidates, AICc ranking, 99 refitted KS replicates, and 10-/200-year plausibility checks. The rainfall 10-year level must lie in 20-800 mm and its 200-year level cannot exceed 1500 mm. Water-level limits are -1 to 5.5 m at 10 years and at most 8 m at 200 years. If no candidate satisfies these criteria, the design value is unavailable.

GBA annual rainfall uses calendar years with at least 274 valid daily values. Mainland units use the local annual-maximum mean multiplied by a distance-limited regional GEV growth curve and the 1.13 daily-to-moving-24-hour factor. Tide years follow the declared annual/July-September coverage rule; partial years whose observed maximum exceeds the median of qualifying years are explicitly identified. Storm-day checks can exclude additional years.

For annual-maximum distribution $F_{\mathrm{AM}}$,

$$x_T=F_{\mathrm{AM}}^{-1}(1-T^{-1}),\qquad T\in\{10,50,100,200\}.$$

The S2 curve belongs to a tide station and is reused for every unit that shares that station. A short station's design curve is transferred from its declared donor using the median difference of concurrent accepted annual maxima:

$$h_T^{(s)}=h_T^{(d)}+\operatorname{median}_{y\in\mathcal Y_{sd}}\{M_y^{(s)}-M_y^{(d)}\}.$$

The transfer adjusts the donor frequency curve by a concurrent annual-maximum offset. Its uncertainty is not quantified.

For eligible GBA mainland units, local annual-maximum mean $\bar M_u$ scales a regional GEV growth curve $g_T$:

$$R_{T,u}^{24\mathrm h}=1.13\,\bar M_u\,g_T.$$

Regional L-CV and L-skewness are weighted by station record length. Only stations within 60 km are eligible. Regional GEV simulations (500 replicates) approximate the Hosking-Wallis $H_1$ heterogeneity statistic. If $H_1\ge2$, the farthest station is removed iteratively. A remaining single station supplies its own growth curve. The 1.13 factor is an empirical daily-to-moving-24-hour conversion; it is an estimate based on daily records.

Alviso's transfer uses accepted concurrent annual maxima from 2019 onward. Other transfer pairs use all concurrent accepted annual samples. Reviewed annual observations can be supplied with explicit selection flags; their files are hashed with the other inputs.

## 04. Fit event margins and copulas

`step04_copula` fits both event margins, ranks dependence models and records bootstrap diagnostics for each unit.

### Marginal distributions

The single-family pool contains Gamma, Lognormal, Weibull, Exponential, GEV, GPD, Log-Pearson III, Log-logistic, Inverse Gaussian, Generalized Gamma, Pearson III, and Burr XII. A two-component lognormal mixture additionally enters the event rainfall candidate pool. The daily configuration does not supply a mixture KS acceptance p value; an eligible single-family fit is preferred. Marginal estimation uses SciPy; mixture estimation uses scikit-learn.

The rainfall event margin models the excess $R-R_{95}$. Its offset restores physical rainfall in CDF, inverse-CDF, and density evaluations. Water level is fitted on the same selected event population.

$$
\operatorname{AICc}=2k-2\ell+\frac{2k(k+1)}{n-k-1}.
$$

Eligibility checks reject nonfinite likelihoods, singular endpoint density modes, inadmissible shape parameters, and weakly identified shifts. Annual candidate parameters exceeding $10^4$ in absolute value are rejected as numerically degenerate. The GBA 10a water-level event margin uses this nondegenerate annual-family pool with event-sample likelihood and KS selection; its dependence family is Gaussian.

### Goodness of fit

Hourly event margins use 49 refitted KS replicates for screening and 999 for final sequential testing in AICc order. Candidates with screening $p<0.02$ do not proceed. GBA daily event margins use 49 refitted Monte Carlo replicates. The mixture KS procedure refits its parameters in every bootstrap sample.

Dependence candidates include Independence, Gaussian, Student, Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7, BB8, and Tawn, with eligible rotations. There are 36 family-rotation combinations. Parameters are estimated by IFM with maximum likelihood in pyvinecopulib. Sequential joint-CDF CvM tests refit both margins and the fixed copula family in every replicate.

$$
S_n=n\,\frac1n\sum_{i=1}^{n}\left\{C_n(\mathbf u_i)-C_{\widehat\theta}(\mathbf u_i)\right\}^2,
\qquad p=\frac{1+\sum_{b=1}^{B_{\mathrm{valid}}}\mathbf1(S_n^{(b)}\ge S_n)}{1+B_{\mathrm{valid}}}.
$$

The joint test uses 99 replicates for San Francisco and Tokyo and 49 for GBA. At least 90% must be valid, and $p\ge0.05$ is required. Formal compound estimates also require both marginal tests to pass, at least 30 event pairs, and at least 10 effective paired years. Otherwise results remain diagnostic. These selection-based tests are not independent external validation.

## 05. Calculate and export designs

`step05_export` calibrates the Kendall contour, selects its physical-density maximum and exports the design table and report. It saves calibration diagnostics, contour points, model status and the analysis settings alongside the figures.

Let $U=F_R(R)$, $V=F_H(H)$, and $K_C(c)=\Pr[C(U,V)\le c]$. For selected-event rate $\lambda$,

$$T_K(c)=\frac{1}{\lambda[1-K_C(c)]},\qquad c_T=K_C^{-1}\!\left(1-\frac{1}{\lambda T}\right).$$

The critical joint-CDF score is calibrated with randomized quasi-Monte Carlo simulation. San Francisco and Tokyo use 1,048,576 calibration samples and an independent equally sized validation sample; relative exceedance-probability error must remain below 5%. GBA uses 262,144 calibration samples with its fixed daily-event numerical settings.

On $C(F_R(r),F_H(h))=c_T$, choose the maximum physical joint-density point:

$$
f_{R,H}(r,h)=c_{RH}(F_R(r),F_H(h))f_R(r)f_H(h).
$$

The hourly workflow combines linear and tail-dense global scans with refinement at all local density peaks. The GBA daily workflow uses its 80-point quantile scan and bounded local refinement. Layer residuals, candidate optima, and calibration metadata are saved. AND and OR exceedance probabilities define different return periods.

### Model scenarios

The scenario matrix has three hazard classes, four return periods, and seven climates: current, SSP1-2.6/SSP2-4.5/SSP3-7.0 at 2060 and 2100. These produce 84 cases. Six C83-based peak-phase cases complete the matrix.

| Hazard | Rainfall boundary | Coastal boundary |
|:--|:--|:--|
| R | S1 rainfall times the prescribed climate multiplier | Typical astronomical tide plus sea-level increment |
| S | No rainfall | S2 peak in model datum plus sea-level increment |
| RS | Kendall rainfall times the prescribed multiplier | Kendall peak in model datum plus sea-level increment |

The regular rainfall/tide peak occurs at hour 13 for San Francisco and hour 12 for GBA and Tokyo. C85-C90 shift the tide peak by -9, -6, -3, +3, +6, and +9 hours relative to the rain peak. The whole coastal water-level process is shifted, preserving its shape.

Climate multipliers, sea-level increments, observed tide templates, rainfall profiles, and river/baseflow boundary records are private inputs. `jrpa cases` creates scenario definitions without those values. `jrpa.scenarios.case_designs` combines supplied climate inputs with locally generated statistical design tables.

### Model vertical reference

| Region/source system | Model conversion |
|:--|:--|
| GBA mainland yearbook | $H_{\mathrm{model}}=H_{1985}-0.744\;\mathrm m$ |
| GBA Hong Kong | $H_{1985}=H_{\mathrm{mCD}}-0.879\;\mathrm m$; $H_{\mathrm{model}}=H_{1985}-0.573\;\mathrm m$ |
| San Francisco Bay | $H_{\mathrm{model}}=H_{\mathrm{NAVD88}}$ |
| Tokyo Bay | $H_{\mathrm{model}}=H_{\mathrm{TP}}$ |

GBA conversions retain the source/model engineering reference ties. They are not an independent geodetic validation of the mainland-to-Hong-Kong offset. Historical water levels retain their trends. Future sea-level increments must carry their stated scenario, epoch, and baseline; they are not a datum conversion.

## Numerical settings

Python 3.12 is the reference interpreter. The package pins NumPy 2.2.6, pandas 2.3.3, SciPy 1.15.3, scikit-learn 1.7.2, Matplotlib 3.11.1, and pyvinecopulib 0.7.5. Seed values, candidate pools, bootstrap counts, and numerical search settings are serialized with every run.

| Setting | GBA | San Francisco / Tokyo |
|:--|--:|--:|
| Event margin KS replicates | 49 | 49 screening; 999 final |
| Joint bootstrap replicates | 49 | 99 |
| Annual margin KS replicates | 99 | 99 |
| Kendall calibration samples | 262,144 | 1,048,576 |
| Independent Kendall validation samples | Not part of daily configuration | 1,048,576 |
| Event fitting seed | 20261011 | 202610052, with unit-specific offsets |
| Annual fitting seed | 20261015 | 20261015 |



The artificial demonstration uses reduced sampling counts and is labelled synthetic. The production settings above apply to study analyses.

## Interpretation

Kendall return periods apply to the selected rainfall-triggered population. S1 and S2 have different annual-maximum populations. A shared tide station has one S2 curve; different rainfall populations can yield different compound points at that station. Separate unit designs do not define a simultaneous bay-wide return period. Historical water levels are not detrended. Parameter-uncertainty intervals are not supplied.

Passing the stated goodness-of-fit and record-length criteria does not establish tail accuracy. Long return periods require extrapolation. Datum uncertainty, short-record transfer, daily-to-24-hour conversion and regional heterogeneity should be considered when interpreting design levels.

## Code locations

| Step | Folder | Supporting modules in that folder |
|--:|:--|:--|
| 01 | `step01_prepare/run.py` | `common`: `io`, `providers`, `datums`, `quality`, `annual_samples` |
| 02 | `step02_events/run.py` | `events`, `event_kernel`, `coverage` |
| 03 | `step03_annual/run.py` | `annual`, `regional`, `annual_fit` |
| 04 | `step04_copula/run.py` | `marginals`, `daily`, `dependence` |
| 05 | `step05_export/run.py` | `kendall`, `reporting`; figure style is `common/plotting` |

`common/` holds settings and checks used by more than one step. It is not an extra step. Study-area maps are already drawn; their drawing program is not part of this release.
