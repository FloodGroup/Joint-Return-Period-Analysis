# Joint Return-Period Analysis

Rainfall and coastal water-level frequency analysis for three bay regions, with annual single-variable designs and rainfall-triggered compound designs at 10, 50, 100 and 200 years.

| Study area | Folder | Units | Observations | Water-level datum |
|:--|:--|--:|:--|:--|
| Greater Bay Area | [gba/](gba/README.md) | 13 | Daily | CN1985 |
| San Francisco Bay | [san_francisco/](san_francisco/README.md) | 11 | Hourly; intact 6 h rainfall supplements | NAVD88 |
| Tokyo Bay | [tokyo/](tokyo/README.md) | 13 | Hourly | Tokyo Peil |

Run the region you need. That region's `run_analysis.py` is the main program:

```text
gba/run_analysis.py
san_francisco/run_analysis.py
tokyo/run_analysis.py
```

Each of those scripts calls the same five folders, in order. Open the script to see the five calls. `shared/jrpa/main.py` is that same sequence in one place.

```text
gba/                     # Greater Bay Area
san_francisco/           # San Francisco Bay
tokyo/                   # Tokyo Bay
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

The same entry files are present in all three regional folders. The five steps are maintained once, under `shared/jrpa/`.

## Run an analysis

Use Python 3.12 and supply a private observation manifest as described in [Inputs and sources](shared/inputs.md).

```sh
python -m pip install -e "."
python tokyo/run_analysis.py --manifest ../private/tokyo/manifest.json --output ../outputs/tokyo
```

Choose the regional folder, review its `settings.py` and `stations.json`, then run its entry script. The five calls are explicit in each script:

| Step | Operation | Folder |
|--:|:--|:--|
| 01 | Prepare observations, quality rules and vertical references | [step01_prepare](shared/jrpa/step01_prepare/run.py) |
| 02 | Merge rainfall sources and select independent event pairs | [step02_events](shared/jrpa/step02_events/run.py) |
| 03 | Fit annual rainfall (S1) and water-level (S2) distributions | [step03_annual](shared/jrpa/step03_annual/run.py) |
| 04 | Fit event margins and copula dependence | [step04_copula](shared/jrpa/step04_copula/run.py) |
| 05 | Calculate Kendall design points and export tables and figures | [step05_export](shared/jrpa/step05_export/run.py) |

The output directory contains `design_values.csv`, diagnostics, source usage and `results.html`. Add `--annual-only` for S1 and S2 alone. Run a small artificial example with `jrpa demo --region tokyo --output ../synthetic-tokyo`.

## Read the documentation

Each regional README includes the drawn maps, station assignments and observation coverage. [Methods](shared/methods.md) gives the equations; [Inputs and sources](shared/inputs.md) describes the observation files and provider references.

To build local HTML pages with typeset equations:

```sh
python -m pip install -e ".[docs]"
python build_docs.py
```

Open `gba/site/index.html`, `san_francisco/site/index.html` or `tokyo/site/index.html`. Use `python build_docs.py --region tokyo` to build one region. The study-area figures are already in each region's `figures/` folder.

## Availability and citation

The repository contains code, station metadata and study-area and coverage figures. Observation records, source GIS, fitted parameters and design-value tables remain private. Obtain observations from the providers listed in the input guide.

Software is licensed under [MIT](LICENSE); third-party notices are in [NOTICE](NOTICE). [CITATION.cff](CITATION.cff) supplies citation metadata. Report issues with the region, processing step, environment and a reproducible artificial example.
