# Analysis workflow

Run the entry script for [GBA](../GBA/run_analysis.py), [SanF](../SanF/run_analysis.py) or [Tokyo](../Tokyo/run_analysis.py). Each calls the five shared steps below, using the region's `settings.py` and `stations.json`. [Methods](methods.md) explains the statistical definitions, equations, regional rules and references. [Inputs and sources](inputs.md) describes the observation files.

| Step | Purpose | Code |
|--:|:--|:--|
| 01 | Prepare observations and vertical references | [step01_prepare/run.py](jrpa/step01_prepare/run.py) |
| 02 | Merge rainfall and construct event samples | [step02_events/run.py](jrpa/step02_events/run.py) |
| 03 | Perform separate univariate rainfall and water-level analyses | [step03_annual/run.py](jrpa/step03_annual/run.py) |
| 04 | Fit event margins and copula dependence | [step04_copula/run.py](jrpa/step04_copula/run.py) |
| 05 | Calculate design combinations and scenario tables | [step05_export/run.py](jrpa/step05_export/run.py) |

## 01. Prepare observations

Read the manifest, validate time axes and units, apply quality exclusions and convert water levels to CN1985 (GBA), NAVD88 (SanF) or Tokyo Peil (Tokyo). Prepare water-level annual samples and load GBA regional rainfall records. Datum offsets carry validity intervals and evidence sources.

Method: [Observation sources and screening](methods.md#01-prepare-observations).

## 02. Construct event samples

Merge rainfall in source-priority order, retaining observed zeros and recording each station's actual contribution. Prepare annual rainfall maxima from their separately screened records. For joint analysis, pair rainfall above its 95th-percentile threshold with the highest water level in the same storm window (SanF and Tokyo) or calendar day (GBA). Decluster the sample and calculate the event rate from effective paired exposure.

Method: [Event windows, thresholds and exposure](methods.md#02-construct-event-samples).

## 03. Perform univariate return-period analysis

Estimate S1 rainfall and S2 coastal water level at 10, 50, 100 and 200 years from independently screened annual maxima. At-site models use 12 candidate families, AICc ranking and a refitted KS test. Mainland GBA rainfall uses regional L-moments and a GEV growth curve. Specified short tide records use donor curves shifted by concurrent annual-maximum differences. Each S2 curve is fitted once per gauge.

Method: [Annual frequency curves and regional analysis](methods.md#03-fit-annual-frequency-curves).

## 04. Fit event margins and copulas

Fit rainfall excess and coincident water level using the single-family pool, with a two-component lognormal mixture also eligible for hourly event rainfall. Fit dependence from 36 copula family–rotation combinations. Rank candidates by AICc, test margins by refitted KS procedures and test the joint model by a bootstrap Cramér–von Mises procedure that refits both margins and the copula.

Method: [Candidate pools, estimation and goodness of fit](methods.md#04-fit-event-margins-and-copulas).

## 05. Calculate and export designs

Calibrate Kendall return-period contours using the fitted copula and event rate. On each contour, select the rainfall–water-level combination with the highest joint density in physical units. Combine these designs and S1/S2 with supplied climate inputs to define rainfall-dominated, storm-surge-dominated and compound scenarios.

Method: [Kendall designs, scenarios and model datums](methods.md#05-calculate-and-export-designs).

[jrpa.main.run](jrpa/main.py) executes this same sequence for `jrpa run` and `jrpa demo`. Common file handling, settings and plotting functions are in `jrpa/common/`. Maps and observation-coverage figures accompany each regional README. The full [reference list](methods.md#references) is maintained with the methods.
