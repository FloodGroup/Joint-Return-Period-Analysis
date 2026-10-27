# Joint Return-Period Analysis

This repository analyzes rainfall and coastal water-level extremes in the Greater Bay Area, San Francisco Bay and Tokyo Bay. It performs separate univariate return-period analyses for rainfall and coastal water level, and a copula-based joint return-period analysis of the two variables. The design scenarios cover rainfall-dominated events, storm-surge-dominated events and compound rainfall鈥搒torm-surge events at 10-, 50-, 100- and 200-year return periods. Compound combinations are selected on Kendall return-period contours using the maximum joint-density criterion.

| Study area | Folder | Units | Observations | Water-level datum |
|:--|:--|--:|:--|:--|
| Greater Bay Area | [GBA/](GBA/README.md) | 13 | Daily | CN1985 |
| San Francisco Bay | [SanF/](SanF/README.md) | 11 | Hourly; intact 6 h rainfall supplements | NAVD88 |
| Tokyo Bay | [Tokyo/](Tokyo/README.md) | 13 | Hourly | Tokyo Peil |

## Analysis units

An analysis unit links a coastal drainage area to a representative rainfall series and an assigned coastal water-level series. Drainage outlets and receiving waters vary around each bay, as do station locations and observation coverage. Units retain these spatial relationships when constructing the three design-event classes.

The partitions follow drainage divides, river outlets and coastal boundary sectors. Several subcatchments can form one unit, while one geographic polygon can have separate analysis units for different rainfall and tide assignments. GBA has 13 analysis units on 11 geographic polygons; SanF has 11 groups of 31 subcatchments; Tokyo has 13 groups of 24 subcatchments. [Methods](shared/methods.md#analysis-units) explains the partitioning principles, and each regional README lists the component areas and tide assignments.

## Repository layout

Each region has a `run_analysis.py` entry point:

```text
GBA/run_analysis.py
SanF/run_analysis.py
Tokyo/run_analysis.py
```

Each entry runs the five steps in `shared/jrpa/`. The command-line interface uses the same sequence in `shared/jrpa/main.py`.

```text
GBA/                     # Greater Bay Area
SanF/                    # San Francisco Bay
Tokyo/                   # Tokyo Bay
    run_analysis.py      # Main program for this region
    settings.py          # Analysis parameters and station rules
    stations.json        # Station IDs, coordinates and source priorities
    README.md            # Regional methods, stations and coverage
    figures/             # Study-area and coverage figures, already drawn
    site/                # Local HTML documentation, generated on demand
shared/
    jrpa/
        main.py          # The five steps, in order
        step01_prepare/  # Observations, quality rules, vertical references
        step02_events/   # Rainfall merge and independent event pairs
        step03_annual/   # Annual rainfall S1 and water-level S2
        step04_copula/   # Event margins and copula dependence
        step05_export/   # Kendall design points, tables and report
        common/          # Settings, files, datums and quality checks
    methods.md           # Equations and processing steps
    inputs.md            # Observation formats and data sources
    assets/              # Documentation styles and equation renderer
```

## Run an analysis

Use Python 3.12 and supply a private observation manifest as described in [Inputs and sources](shared/inputs.md).

```sh
python -m pip install -e "."
python Tokyo/run_analysis.py --manifest ../private/tokyo/manifest.json --output ../outputs/tokyo
```

Review the region's `settings.py` and `stations.json` before running its entry script. The five steps are:

| Step | Operation | Folder |
|--:|:--|:--|
| 01 | Prepare observations, quality rules and vertical references | [step01_prepare](shared/jrpa/step01_prepare/run.py) |
| 02 | Merge rainfall sources and select independent event pairs | [step02_events](shared/jrpa/step02_events/run.py) |
| 03 | Fit annual rainfall (S1) and water-level (S2) distributions | [step03_annual](shared/jrpa/step03_annual/run.py) |
| 04 | Fit event margins and copula dependence | [step04_copula](shared/jrpa/step04_copula/run.py) |
| 05 | Calculate Kendall design points and export tables and figures | [step05_export](shared/jrpa/step05_export/run.py) |

The output directory contains `design_values.csv`, diagnostics, source usage and `results.html`. Add `--annual-only` for S1 and S2 alone. Run a small artificial example with `jrpa demo --region tokyo --output ../synthetic-tokyo`.

## Read the documentation

Each regional README contains maps, station assignments and observation coverage. [Methods](shared/methods.md) gives the equations; [Inputs and sources](shared/inputs.md) describes the observation files and their sources.

To build local HTML pages with typeset equations:

```sh
python -m pip install -e ".[docs]"
python build_docs.py
```

Open `GBA/site/index.html`, `SanF/site/index.html` or `Tokyo/site/index.html`. Use `python build_docs.py --region tokyo` to build one region. The command-line region identifiers are `gba`, `san_francisco` and `tokyo`.

## Availability and citation

The repository contains code, station metadata and study-area and coverage figures. Observation records, source GIS, fitted parameters and design-value tables remain private. Obtain observations from the providers listed in [Inputs and sources](shared/inputs.md).

Software is released under the MIT License ([LICENSE](LICENSE)). Local HTML documentation bundles KaTeX 0.16.22 (MIT) under `shared/assets/vendor/katex/`.

Please cite this software as:

> Yang, R., Long, Z., Jiang, L., & Duan, H.-F. (2026). *Joint Return-Period Analysis* (Version 1.0.0) [Computer software]. The State Key Laboratory of Climate Resilience for Coastal Cities (SKL-CRCC), The Hong Kong Polytechnic University. https://github.com/FloodGroup/Joint-Return-Period-Analysis

```bibtex
@software{Yang_Joint_Return_Period_Analysis_2026,
  author       = {Yang, Ruiyi and Long, Zhiyong and Jiang, Long and Duan, Huan-Feng},
  title        = {Joint Return-Period Analysis},
  version      = {1.0.0},
  year         = {2026},
  organization = {The State Key Laboratory of Climate Resilience for Coastal Cities (SKL-CRCC), The Hong Kong Polytechnic University},
  url          = {https://github.com/FloodGroup/Joint-Return-Period-Analysis}
}
```

Report issues with the region, processing step, environment and a reproducible artificial example.
