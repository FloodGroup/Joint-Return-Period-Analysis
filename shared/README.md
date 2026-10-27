# Analysis workflow

Run the entry script for the study area:

```text
GBA/run_analysis.py
SanF/run_analysis.py
Tokyo/run_analysis.py
```

Each script calls the five steps in `jrpa/`:

| Step | Folder | Operation | Entry point |
|--:|:--|:--|:--|
| 01 | `jrpa/step01_prepare/` | Read observations, apply quality controls and convert water-level datums | `run.py` |
| 02 | `jrpa/step02_events/` | Merge rainfall and select independent rainfall and water-level events | `run.py` |
| 03 | `jrpa/step03_annual/` | Fit annual rainfall S1 and water-level S2 distributions | `run.py` |
| 04 | `jrpa/step04_copula/` | Fit event margins and copula dependence | `run.py` |
| 05 | `jrpa/step05_export/` | Calculate Kendall design values and export tables and reports | `run.py` |

`jrpa/main.py` runs the same steps for the `jrpa run` and `jrpa demo` commands. Shared settings, file handling, datum conversions and quality controls are in `jrpa/common/`.

Study-area maps are in each region's `figures/` folder. The program used to draw them is not included in this release.
