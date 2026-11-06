# Methods

The analysis defines rainfall-dominated, storm-surge-dominated and compound rainfall–storm-surge design events for the Greater Bay Area (GBA), San Francisco Bay (SanF) and Tokyo Bay (Tokyo). Separate annual-maximum analyses provide the univariate design rainfall (S1) and coastal water level (S2). A bivariate copula describes storm rainfall and the highest water level in the same event window; its Kendall return-period contours provide compound design combinations. All three analyses use target return periods of 10, 50, 100 and 200 years.

Rainfall $R$ is a depth in millimetres, and coastal water level $H$ is an elevation in metres. The marine variable includes astronomical tide and meteorological surge. The methods follow the frequency-analysis description in *Reclaimed land more prone to compound flooding in three Pacific bay areas* and its Supplementary Note 4; scenario definitions follow Supplementary Notes 4–6.

[Units](#analysis-units) · [01 Observations](#01-prepare-observations) · [02 Events](#02-construct-event-samples) · [03 Univariate analysis](#03-fit-annual-frequency-curves) · [04 Copula analysis](#04-fit-event-margins-and-copulas) · [05 Designs](#05-calculate-and-export-designs) · [References](#references)

## Analysis units

An analysis unit links a coastal drainage area to a representative rainfall record and the tide gauge serving its receiving channel or coastal sector. This spatial assignment allows rainfall designs to vary between drainage systems while maintaining a consistent coastal boundary. Several units can share a tide gauge; each retains its own rainfall sample and compound-event model, while the shared gauge has one S2 frequency curve.

Drainage divides, river outlets and coastal boundary sectors determine the geographic partitions. The terrain workflows use Priority-Flood depression filling [1] and D8 routing, with river-network constraints in flat terrain. Small coastal areas are grouped by drainage destination and land adjacency. These polygons are supplied by the model setup; the frequency-analysis steps use their fixed station assignments.

| Region | Geographic basis | Statistical units |
|:--|:--|:--|
| GBA | 30 m GLO-30 terrain, GEOGloWS v2/TDX-Hydro river constraints and D8 drainage destinations; 11 outlet-based coastal polygons | 13 units; polygons 3 and 10 each have two rainfall–tide assignments |
| SanF | 30 m USGS 3DEP terrain, NHDPlusV2 channels and WBD HUC12 constraints; 31 coastal subcatchments | 11 coastal-drainage groups |
| Tokyo | 30 m terrain aggregated from GSI elevation data, MLIT W05 rivers and W07 basin boundaries; 24 coastal subcatchments | 13 coastal-drainage groups |

GBA units 03a/03b and 10a/10b share geographic polygons but use different observation assignments. In `stations.json`, `unit` identifies the statistical assignment and `geometry_id` the mapped polygon. The regional documentation lists component subcatchments, station IDs, source priorities and actual observation coverage: [GBA](../GBA/README.md), [SanF](../SanF/README.md) and [Tokyo](../Tokyo/README.md).

## 01. Prepare observations

### Observation sources and time resolution

| Region | Rainfall | Coastal water level | Analysis clock |
|:--|:--|:--|:--|
| GBA | CMA automatic gauges, aggregated by calendar day; HKO and GEO daily gauges; long-record daily gauges for regional frequency analysis | Chinese hydrological yearbooks and HKO records of daily maximum water level | Local calendar days, Hong Kong time |
| SanF | GHCNh, NCEI-COOP and CDEC gauges; hourly observations supplemented by intact six-hour accumulation reports | NOAA and USGS gauges | Hourly; fixed PST for annual grouping |
| Tokyo | JMA gauges | Japanese gauges in GESLA-3 and JODC records | Hourly; JST for annual grouping |

Rainfall sources are merged in the declared distance-priority order. A valid observation is retained, including zero rainfall; a missing value is filled from the next available source without rescaling. In the hourly workflow, observed hourly values take priority over six-hour reports. A six-hour total is accepted only as an intact block when none of its hours already contains an accepted observation. This avoids counting overlapping accumulations twice.

Quality exclusions are applied before sample construction. Annual rainfall exclusions are applied to the annual-analysis merge independently of the event-analysis merge. Each accepted accumulation retains its source station, allowing the station tables and coverage figures to describe the data actually used. File formats, source links and manifest fields are documented in [Inputs and sources](inputs.md).

### Vertical reference and record screening

Water levels are converted to CN1985 in GBA, NAVD88 in SanF and Tokyo Peil (TP) in Tokyo before statistical fitting. Each datum conversion has an offset, a validity interval and an evidence source. Epochs with an unknown offset remain missing. The GBA model-level conversion is applied after fitting, as described under [Model vertical reference](#model-vertical-reference).

Documented sensor errors and missing or distorted storm peaks are screened explicitly. Historical water-level trends are retained. Annual rainfall and water-level records are screened separately, so their univariate analyses can use different observation periods. Simultaneous coverage is required for the compound-event sample.

## 02. Construct event samples

The joint sample represents heavy-rainfall events and their coincident coastal water levels. Rainfall must exceed its selection threshold; water level has no additional exceedance threshold. This retains both high and moderate coastal levels during heavy rain, as specified in Supplementary Note 4.4. Water-level-driven events are represented separately by S2.

### SanF and Tokyo: storm-window pairing

1. Identify positive rainfall accumulations separated by intervals of at least 12 hours. Record whether the intervening dry interval is confirmed by observations; missing rainfall is never converted to zero.
2. For an event shorter than 24 hours, use its complete rainfall total. For a longer event, select its largest complete moving 24-hour total. Let this total be $R_i$ and its selected window be $W_i$.
3. Pair $R_i$ with the highest observed water level in $W_i$, denoted $H_i$. At least 75% of the window's water-level hours must be observed.
4. Compute $R_{95}$ from all rainfall-complete event totals, before rejecting windows for water-level coverage. Retain strict exceedances, $R_i>R_{95}$.
5. Group events whose adjacent water-level peaks are less than 72 hours apart. Retain the largest rainfall in each cluster; an earlier selected window resolves a tie.

The 12-hour separation identifies rainfall events; the subsequent 72-hour clustering controls repeated sampling of the same storm system.

### GBA: daily pairing

The rainfall threshold is the 95th percentile of positive daily rainfall within the declared analysis clock. Days exceeding this threshold are paired with the highest water level of the same calendar day. Adjacent exceedance dates less than three days apart form a cluster. The day with the largest rainfall is retained, with water level resolving rainfall ties. Explicitly invalid selected pairs are removed without selecting a replacement day from that cluster.

### Effective exposure and event rate

The event rate converts a probability within the selected storm population to a recurrence interval in years. Let $N_{\mathrm{paired}}$ be the number of hourly or daily intervals with simultaneous valid rainfall and water-level observations, and $n$ the number of retained events. Then

```math
Y_{\mathrm{eff}}=\frac{N_{\mathrm{paired}}}{365.2425\,s},
\qquad \lambda=\frac{n}{Y_{\mathrm{eff}}},
```

where $s=24$ for hourly data and $s=1$ for daily data. Thus $Y_{\mathrm{eff}}$ is effective paired exposure in years and $\lambda$ is the mean number of retained events per year. Missing intervals reduce exposure even when they fall between the first and last observation dates.

## 03. Fit annual frequency curves

### Annual samples and univariate return periods

S1 and S2 are separate univariate return-period analyses. S1 uses the annual maximum rainfall of each unit's full composite record. S2 uses the annual maximum water level of each gauge's full record. Neither sample is restricted to the years available for joint analysis.

| Record | Annual grouping and acceptance |
|:--|:--|
| SanF and Tokyo rainfall | Complete moving 24-hour totals; October–September water years in fixed PST and JST, respectively; at least 75% valid hours |
| SanF and Tokyo water level | Highest hourly level in each October–September water year; at least 75% valid hours |
| GBA rainfall | Calendar-year maximum daily rainfall; at least 274 valid days |
| GBA water level | Calendar-year maximum; at least 75% annual coverage or 80% July–September coverage. Additional partial years are retained when their observed maximum reaches the median of qualifying annual maxima; these years are recorded explicitly. Storm-day exclusions are then applied. |

For an annual-maximum cumulative distribution $F_{\mathrm{AM}}$, a return period of $T$ years corresponds to annual exceedance probability $1/T$:

```math
x_T=F_{\mathrm{AM}}^{-1}(1-T^{-1}),
\qquad T\in\{10,50,100,200\}.
```

Here $x_T$ is rainfall depth for S1 or water-level elevation for S2. At-site fits use maximum likelihood and the 12 single-family candidates listed in [Marginal distributions](#marginal-distributions). Candidates are ranked by AICc [2] and tested in that order. The first candidate passing the refitted KS test ($p\ge0.05$), numerical checks and the following return-level checks is retained.

| Variable | 10-year value | 200-year value |
|:--|:--|:--|
| Rainfall | 20–800 mm | At most 1,500 mm |
| Water level, in analysis datum | −1 to 5.5 m | At most 8 m |

### GBA regional rainfall frequency analysis

Mainland GBA units 01–08, including 03a and 03b, use an index-rainfall analysis based on regional L-moments [3]. Hong Kong units 09, 10a, 10b and 11 use at-site fits. The regional procedure separates the local rainfall magnitude from the dimensionless frequency-curve shape:

```math
R_{T,u}^{24\mathrm h}=a\,\bar M_u\,g_T,
\qquad a=1.13.
```

Here $\bar M_u$ is the mean accepted annual maximum at unit $u$, and $g_T$ is a dimensionless GEV growth factor. Eligible donor gauges lie within 60 km. Their L-CV and L-skewness are averaged with weights proportional to record length to fit the regional GEV. Heterogeneity is assessed from the dispersion of donor L-CVs using 500 simulations from the fitted regional GEV, a study implementation of the Hosking–Wallis $H_1$ approach. If $H_1\ge2$, the farthest donor is removed and the assessment repeated. A remaining single donor supplies its own growth curve.

The factor 1.13 converts fixed daily rainfall to estimated moving 24-hour rainfall [4, 5]. Supplementary Note 4.3 supports this choice with a median moving-to-fixed ratio of 1.125 from 118 annual-maximum windows at 34 hourly gauges. The same factor is applied to GBA at-site S1 values and, after joint design-point selection, to compound rainfall. It is applied once; SanF and Tokyo use their observed moving totals directly.

### Short-record water-level transfer

A specified short gauge uses the donor's S2 curve shifted by the median difference of concurrent accepted annual maxima:

```math
h_T^{(s)}=h_T^{(d)}+
\mathrm{median}_{y\in\mathcal Y_{sd}}
\{M_y^{(s)}-M_y^{(d)}\}.
```

Superscripts $s$ and $d$ denote the short and donor gauges; $\mathcal Y_{sd}$ contains their common accepted years. This preserves the donor curve's shape while adjusting its elevation to the short gauge.

| Region | Short gauge | Donor |
|:--|:--|:--|
| GBA | Dajiuzhou; Neilingding | Sanzao; Chiwan, respectively |
| SanF | Alviso | Redwood City; accepted common years from 2019 onward |
| Tokyo | Kawasaki; Kurihama | Tokyo-Harumi; Aburatsubo, respectively |

Other transfers use all concurrent accepted annual samples. Each gauge's resulting S2 curve is assigned to every unit served by that gauge.

## 04. Fit event margins and copulas

### Marginal distributions

Copula analysis separates the marginal distributions of event rainfall and coastal water level from their dependence [6]. The event rainfall margin is fitted to the positive excess $E=R-R_{95}$. The threshold is restored in all evaluations of the physical rainfall CDF, inverse CDF and density. The water-level margin is fitted to $H$ from the same selected events. These event distributions are separate from the annual-maximum distributions used for S1 and S2.

| Candidate pool | Families |
|:--|:--|
| Single-family pool | Gamma, lognormal, Weibull, exponential, generalized extreme value, generalized Pareto, log-Pearson III, log-logistic, inverse Gaussian, generalized gamma, Pearson III and Burr XII |
| Additional hourly event-rainfall candidate | Two-component lognormal mixture, fitted by expectation–maximization |

Single-family parameters are estimated by maximum likelihood. GBA uses the single-family pool for both event margins. Eligibility checks reject nonfinite likelihoods, singular endpoint density modes, inadmissible shapes and weakly identified shifts. Annual candidate parameters exceeding $10^4$ in absolute value are rejected as numerically degenerate. GBA unit 10a applies the annual-family eligibility rules to its event water-level sample and uses a Gaussian copula.

Eligible candidates are ranked by the corrected Akaike information criterion [2]:

```math
\mathrm{AICc}=2k-2\ell+\frac{2k(k+1)}{n-k-1},
```

where $\ell$ is the maximized log-likelihood, $k$ the number of fitted parameters and $n$ the sample size. The same ranking criterion is used for the annual fits and copula candidates, with the likelihood and parameter count of the respective model.

### Copula families and estimation

For fitted marginal CDFs $F_R$ and $F_H$, the joint CDF is

```math
F_{R,H}(r,h)=C_\theta\{F_R(r),F_H(h)\},
\qquad u_i=F_R(R_i),\quad v_i=F_H(H_i).
```

The inference-functions-for-margins procedure first estimates both margins, then estimates the copula parameters $\theta$ by maximum likelihood on the transformed pairs [6]. Copula fitting uses pyvinecopulib [9].

| Dependence family | Rotations |
|:--|:--|
| Independence, Gaussian, Student-t, Frank | 0° |
| Clayton, Gumbel, Joe, BB1, BB6, BB7, BB8, Tawn | 0°, 90°, 180°, 270° |

The pool therefore contains 36 family–rotation combinations. Rotations admit different directions of tail dependence, including negative dependence. The independence model is an explicit candidate. Copulas are tested in AICc order, and the first passing the joint goodness-of-fit test is selected.

### Goodness of fit

The marginal Kolmogorov–Smirnov statistic measures the largest difference between the empirical CDF $F_n$ and the fitted CDF:

```math
D_n=\sup_x\left|F_n(x)-F_{\widehat\psi}(x)\right|.
```

Its null distribution is estimated by parametric simulation with parameters refitted to every replicate [8]. This accounts for estimating $\psi$ from the tested sample. A lognormal mixture is also fully refitted in every replicate. Annual fits use 99 replicates. Hourly event fits use 49 screening replicates; candidates with screening $p\ge0.02$ proceed, in AICc order, to 999-replicate testing. GBA event fits use 49 replicates. Final marginal acceptance requires $p\ge0.05$.

The joint test uses a Cramér–von Mises discrepancy evaluated at the fitted marginal transforms. Define the empirical joint CDF of these transforms as

```math
C_n(u,v)=\frac1n\sum_{j=1}^{n}
\mathbf1\{u_j\le u,\ v_j\le v\}.
```

Then the statistic is

```math
S_n=\sum_{i=1}^{n}
\left[C_n(u_i,v_i)-C_{\widehat\theta}(u_i,v_i)\right]^2.
```

The bootstrap follows the empirical-CDF testing principle reviewed by Genest et al. [7]. In this implementation, each replicate simulates paired observations from the fitted joint model, refits both marginal distributions, transforms the simulated sample with those refitted margins, and refits the same copula family and rotation. It then recomputes $S_n$. Family selection is held fixed within a replicate. The Monte Carlo p value is

```math
p=\frac{1+\sum_{b=1}^{B_{\mathrm{valid}}}
\mathbf1\{S_n^{(b)}\ge S_n\}}{1+B_{\mathrm{valid}}},
```

where $B_{\mathrm{valid}}$ counts successful replicates. The joint test uses 49 replicates for GBA and 99 for SanF and Tokyo; at least 90% must be valid and $p\ge0.05$. The compound-fit acceptance flag also requires both marginal tests to pass, at least 30 event pairs and 10 effective paired years. The run records these checks separately from the fitted model.

## 05. Calculate and export designs

### Kendall return period

The Kendall construction assigns a recurrence interval to a critical level of the joint CDF [10]. Let $(U,V)$ follow the fitted copula and define

```math
K_C(c)=\Pr\{C(U,V)\le c\}.
```

The probability of a selected event lying above critical score $c$ is $1-K_C(c)$. Multiplying this by the event rate gives its mean annual occurrence rate. Consequently,

```math
T_K(c)=\frac{1}{\lambda[1-K_C(c)]},
\qquad c_T=K_C^{-1}\!\left(1-\frac{1}{\lambda T}\right).
```

All pairs on $C(F_R(r),F_H(h))=c_T$ share this Kendall return period. The contour is calibrated by simulating from the fitted copula, evaluating its CDF scores and taking the score quantile at $1-1/(\lambda T)$. GBA uses 262,144 quasi-random samples. SanF and Tokyo use 1,048,576 calibration samples and an equally sized independent validation sample; the relative error in the validation exceedance probability must be below 5%.

### Representative compound design point

A contour contains many rainfall–water-level combinations. The study selects the point with the highest joint density in physical coordinates [10]:

```math
f_{R,H}(r,h)=c_\theta(F_R(r),F_H(h))\,f_R(r)\,f_H(h),
```

```math
(r_T,h_T)=\underset{C(F_R(r),F_H(h))=c_T}{\mathrm{arg\,max}}
\ f_{R,H}(r,h).
```

Here $c_\theta$ is the copula density, while $f_R$ and $f_H$ are the marginal densities. Including both marginal densities makes the optimization apply to rainfall depth and water-level elevation. For each trial rainfall quantile, the corresponding water-level quantile is solved on the contour. SanF and Tokyo combine linear and tail-dense scans with local refinement at every detected density maximum. GBA uses an 80-point quantile scan and bounded refinement around its best point. The fitted contour residual must not exceed $10^{-8}$. GBA design rainfall is then multiplied by 1.13.

For reference, joint AND and OR recurrence intervals at an event-scale pair $(r,h)$ are

```math
T_{\mathrm{AND}}(r,h)=
\frac{1}{\lambda[1-u-v+C(u,v)]},
\qquad T_{\mathrm{OR}}(r,h)=
\frac{1}{\lambda[1-C(u,v)]},
```

where $u=F_R(r)$ and $v=F_H(h)$. These describe different exceedance sets from the Kendall construction. The RS designs use $T_K$; their individual rainfall and water-level components receive their annual-maximum recurrence intervals from the separate S1 and S2 distributions. Thus an RS100 design is a point on the 100-year Kendall contour, rather than a pairing of the two univariate 100-year levels.

### Model scenarios

The scenario matrix contains three event classes, four return periods and seven climate settings: present conditions and SSP1-2.6, SSP2-4.5 and SSP3-7.0 in 2060 and 2100. This gives 84 cases, C01–C84. Six timing variants of C83, the 100-year compound event under SSP3-7.0 in 2100, complete the 90-case matrix.

| Scenario | Rainfall boundary | Coastal boundary |
|:--|:--|:--|
| R: rainfall-dominated | S1 rainfall multiplied by the climate factor | Typical astronomical tide plus sea-level increment |
| S: storm-surge-dominated | No rainfall | S2 peak in model datum plus sea-level increment |
| RS: compound rainfall–storm-surge | Kendall rainfall multiplied by the climate factor | Kendall water-level peak in model datum plus sea-level increment |

The nominal event window is 24 hours. The regular rainfall and storm-tide peaks coincide at hour 13 in SanF and hour 12 in GBA and Tokyo. C85–C90 shift the coastal process by −9, −6, −3, +3, +6 and +9 hours relative to the rainfall peak, preserving its shape. Where several units share a gauge, the paper's hydrodynamic boundary setup uses the largest of their RS design water levels.

Future scenarios scale rainfall and add the supplied sea-level increment while retaining the present-climate event rate, margins and copula. Rainfall factors, sea-level increments and their baselines are supplied as inputs. `jrpa cases` creates the scenario definitions; `common.scenarios.case_designs` applies the numerical climate inputs to the unit-level statistical designs. The rainfall profiles, coastal hydrographs and spatial boundary assembly described in Supplementary Note 5 belong to the subsequent model-input preparation.

### Model vertical reference

The fitting datum and model datum are retained separately in `H_m` and `H_model_m`. The study uses the following conversions:

| Region/source system | Model conversion |
|:--|:--|
| GBA mainland yearbook | $H_{\mathrm{model}}=H_{1985}-0.744\quad\mathrm{m}$ |
| GBA Hong Kong | $H_{1985}=H_{\mathrm{mCD}}-0.879\quad\mathrm{m}$; $H_{\mathrm{model}}=H_{1985}-0.573\quad\mathrm{m}$ |
| SanF | $H_{\mathrm{model}}=H_{\mathrm{NAVD88}}$ |
| Tokyo | $H_{\mathrm{model}}=H_{\mathrm{TP}}$ |

The two Hong Kong offsets give $H_{\mathrm{model}}=H_{\mathrm{mCD}}-1.452\quad\mathrm{m}$, consistent with Supplementary Note 4.1. GBA model elevations use the study's source-specific mean-sea-level reference ties. Future sea-level increments are added after this conversion.

## Numerical settings

| Setting | GBA | SanF / Tokyo |
|:--|--:|--:|
| Annual-margin KS replicates | 99 | 99 |
| Event-margin KS replicates | 49 | 49 screening; 999 final |
| Joint bootstrap replicates | 49 | 99 |
| Kendall calibration samples | 262,144 | 1,048,576 |
| Independent Kendall validation samples | 0 | 1,048,576 |
| Event fitting seed | 20261011 | 202610052, with unit-specific offsets |
| Annual fitting seed | 20261015 | 20261015 |

Regional settings store these values and the sampling rules. The production environment uses Python 3.12 and the versions pinned in [pyproject.toml](../pyproject.toml), including SciPy 1.15.3 and pyvinecopulib 0.7.5. Run metadata retains the settings and input hashes. The synthetic demonstration uses reduced simulation counts.

## Code locations

Each regional `run_analysis.py` and `jrpa.main.run` calls the same five steps. Regional `settings.py` files specify parameters and station rules; `stations.json` records the spatial and observation assignments.

| Step | Entry point | Main supporting modules |
|--:|:--|:--|
| 01 | [Prepare observations](jrpa/step01_prepare/run.py) | `common/io`, `providers`, `datums`, `quality`, `annual_samples` |
| 02 | [Construct samples](jrpa/step02_events/run.py) | `events`, `event_kernel`, `coverage` |
| 03 | [Univariate analysis](jrpa/step03_annual/run.py) | `annual`, `regional`, `annual_fit` |
| 04 | [Copula analysis](jrpa/step04_copula/run.py) | `marginals`, `daily`, `dependence` |
| 05 | [Designs and reporting](jrpa/step05_export/run.py) | `kendall`, `reporting`; shared styles in `common/plotting` |

## References

1. Barnes, R., Lehman, C. & Mulla, D. (2014). Priority-Flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. *Computers & Geosciences*, **62**, 117–127. [doi:10.1016/j.cageo.2013.04.024](https://doi.org/10.1016/j.cageo.2013.04.024).
2. Hurvich, C. M. & Tsai, C.-L. (1989). Regression and time series model selection in small samples. *Biometrika*, **76**, 297–307. [doi:10.1093/biomet/76.2.297](https://doi.org/10.1093/biomet/76.2.297).
3. Hosking, J. R. M. & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press. [doi:10.1017/CBO9780511529443](https://doi.org/10.1017/CBO9780511529443).
4. Hershfield, D. M. (1961). *Rainfall Frequency Atlas of the United States for Durations from 30 Minutes to 24 Hours and Return Periods from 1 to 100 Years*. Technical Paper No. 40, US Weather Bureau. [NOAA archive](https://www.weather.gov/media/owp/oh/hdsc/docs/TechnicalPaper_No40.pdf).
5. World Meteorological Organization (2009). *Manual on Estimation of Probable Maximum Precipitation (PMP)*. WMO-No. 1045, 3rd ed. [WMO publication catalog](https://community.wmo.int/site/knowledge-hub/programmes-and-initiatives/water-resources-assessment/hydrology-publications).
6. Salvadori, G., De Michele, C., Kottegoda, N. T. & Rosso, R. (2007). *Extremes in Nature: An Approach Using Copulas*. Springer. [doi:10.1007/1-4020-4415-1](https://doi.org/10.1007/1-4020-4415-1).
7. Genest, C., Rémillard, B. & Beaudoin, D. (2009). Goodness-of-fit tests for copulas: A review and a power study. *Insurance: Mathematics and Economics*, **44**, 199–213. [doi:10.1016/j.insmatheco.2007.10.005](https://doi.org/10.1016/j.insmatheco.2007.10.005).
8. SciPy developers. *scipy.stats.goodness_of_fit*, SciPy 1.15.3. [Refitted Monte Carlo goodness-of-fit procedure](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.goodness_of_fit.html).
9. Nagler, T. & Vatter, T. *pyvinecopulib*. [Software archive, doi:10.5281/zenodo.4288292](https://doi.org/10.5281/zenodo.4288292); [version 0.7.5 used by this repository](https://github.com/vinecopulib/pyvinecopulib/tree/v0.7.5).
10. Salvadori, G., De Michele, C. & Durante, F. (2011). On the return period and design in a multivariate framework. *Hydrology and Earth System Sciences*, **15**, 3293–3305. [doi:10.5194/hess-15-3293-2011](https://doi.org/10.5194/hess-15-3293-2011).
