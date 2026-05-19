# Compound Flood Return-Period Analysis

Open-source scripts for estimating **bivariate rainfall–sea level return periods** in three megacity bay areas:

- **GBA** — Greater Bay Area estuary (mainland China + Hong Kong)
- **San Francisco** — San Francisco Bay
- **Tokyo** — Tokyo Bay

---

## 1. Method

Each regional script follows the same steps:

1. **Pair stations** — Assign one or more rainfall gauges to each tide gauge (fixed pairs on the mainland; priority-ordered lists elsewhere; San Francisco uses the six nearest gauges by distance).
2. **Merge and align** — Build a single daily rainfall series from paired gauges (first available date wins by priority); inner-join with daily maximum sea level on calendar date.
3. **Datum correction** — Apply piecewise vertical offsets where tide records use different reference datums (mainland GBA, Tokyo) or a uniform shift (Hong Kong: −0.868 m to 1985 Chinese Height Datum; San Francisco: per-station `datum_offset`).
4. **Define compound events** — Flag days when 24-hour cumulative rainfall and daily maximum sea level both exceed their **95th-percentile** thresholds in the merged series.
5. **De-cluster** — Treat events within **3 days** of each other as one event.
6. **Fit models** — Select marginal distributions (AIC among Gamma, GEV, GP, logistic, normal, Gumbel, Nakagami) and copulas (AIC among Gaussian, Student-t, Clayton, Gumbel, Frank, Joe); estimate by MLE; assess marginals with the Kolmogorov–Smirnov test.
7. **Return periods** — For joint exceedance at levels (r, s):

   \[
   T(R>r,\ S>s)=\frac{1}{1-F_R(r)-F_S(s)+C(F_R(r),F_S(s))}
   \]

   where \(F_R\), \(F_S\) are marginal CDFs and \(C\) is the fitted copula.
8. **Output** — Return-period contour plots; CSV tables of contour axis intersections and maximum-density points.

---

## 2. Data pairing and coverage

Coverage timelines are in [docs/figures/coverage_GBA.png](docs/figures/coverage_GBA.png), [docs/figures/coverage_San_Francisco.png](docs/figures/coverage_San_Francisco.png), and [docs/figures/coverage_Tokyo.png](docs/figures/coverage_Tokyo.png). Blue = tide observations; green = rainfall; red hatching = dates that rainfall actually supplies after priority merge and tide inner join. Percentages are each rainfall gauge's share of the analysable sample.

**Table columns (all regions):**

| Column | Meaning |
|---|---|
| Tide gauge | Tide observation station |
| Rainfall pairing | Gauges used, in priority order (share of analysis sample in parentheses) |
| Tide observations | Span and day count of tide record |
| Analysis sample | Span and day count after tide–rainfall inner join (used in copula fitting) |

### 2.1 GBA estuary — mainland

| Tide gauge | Rainfall pairing | Tide observations | Analysis sample |
|---|---|---|---|
| Chiwan | Tiegang (100%) | 1965–2022 · 13,878 d | 1975–2022 · 10,226 d |
| Sanzao | Sanzao (100%) | 1965–2022 · 15,552 d | 1965–2022 · 15,552 d |
| Nei Lingding | Changjiang (100%) | 2010–2022 · 4,016 d | 2010–2022 · 4,016 d |

Mainland rainfall workbooks treat blank cells as 0 mm daily rainfall.

### 2.2 GBA estuary — Hong Kong

| Tide gauge | Rainfall pairing | Tide observations | Analysis sample |
|---|---|---|---|
| Quarry Bay | Quarry Bay (27.5%) → Kai Tak (0.3%) → Shau Kei Wan (0.1%) → Hong Kong Observatory (72.1%) | 1960–2023 · 22,487 d | 1960–2023 · 22,487 d |
| Tsim Bei Tsui | Tsim Bei Tsui (94.6%) → Wetland Park (1.0%) → Lau Fau Shan (4.4%) | 1983–2023 · 12,847 d | 1985–2023 · 12,281 d |
| Shek Pik | Ngong Ping Fresh Water Reservoir (100%) | 1998–2023 · 9,073 d | 2006–2023 · 6,147 d |
| Tai Mo Wan | Sai Kung (Hong Kong Adventist College) (98.2%) → Shau Kei Wan (1.8%) | 1996–2023 · 9,076 d | 1996–2023 · 8,904 d |
| Tai Po Kau | Tai Po Wong Shiu Chi Secondary School (100%) | 1981–2023 · 15,248 d | 1985–2023 · 13,091 d |
| Waglan Island | Waglan Island (100%) | 1982–2023 · 8,980 d | 1989–2023 · 7,869 d |

### 2.3 San Francisco Bay

Six nearest rainfall gauges are assigned per tide station; gauges contributing less than 2% of the merged sample are omitted from analysis. All tide records use a common vertical datum via station-specific offsets.

| Tide gauge | Rainfall pairing | Tide observations | Analysis sample |
|---|---|---|---|
| San Francisco | San Francisco Downtown (83.0%) → Kentfield (10.3%) → Berkeley (6.6%) | 1898–2021 · 44,331 d | 1898–2021 · 43,259 d |
| Redwood City | Palo Alto (76.5%) → Newark (13.6%) → Woodside Fire Station (7.1%) → Fremont (2.7%) | 1974–2021 · 9,569 d | 1974–2021 · 9,569 d |
| Alameda | Oakland Museum (92.0%) → Oakland Metro Intl AP (5.3%) → San Francisco Downtown (2.6%) | 1976–2021 · 15,366 d | 1976–2021 · 15,366 d |
| Richmond | Richmond (90.8%) → Kentfield (9.2%) | 1996–2021 · 9,065 d | 1996–2021 · 9,065 d |
| Port Chicago | Concord Buchanan Field (53.5%) → Martinez Water Plant (46.3%) | 1979–2021 · 15,038 d | 1979–2021 · 15,038 d |

### 2.4 Tokyo Bay

| Tide gauge | Rainfall pairing | Tide observations | Analysis sample |
|---|---|---|---|
| Chiba | Yokoshibahikari (100%) | 1964–2019 · 20,067 d | 1966–2019 · 19,603 d |
| Tokyo | Tokyo (100%) | 1960–2019 · 20,244 d | 1964–2019 · 19,151 d |
| Yokosuka | Miura (99.9%) → Enoshima (0.1%) | 1960–2019 · 21,565 d | 1976–2019 · 15,976 d |
| Aburatsubo | Miura (100%) | 1932–2019 · 30,430 d | 1976–2019 · 15,875 d |
| Mera | Sakahata (100%) | 1964–2019 · 19,198 d | 1968–2019 · 18,688 d |
| Kawasaki | Haneda (72.4%) → Yokohama (27.6%) | 1967–1996 · 8,771 d | 1967–1996 · 8,771 d |

### 2.5 Data sources

| Region | Tide data | Rainfall data |
|---|---|---|
| GBA mainland | _[TODO: provider URL]_ | _[TODO: provider URL]_ |
| Hong Kong | _[TODO: HKO / CEDD URL]_ | _[TODO: HKO rainfall URL]_ |
| San Francisco | [NOAA CO-OPS](https://tidesandcurrents.noaa.gov/) | [NOAA/NCEI](https://www.ncei.noaa.gov/) |
| Tokyo | _[TODO: JMA / MLIT URL]_ | _[TODO: JMA AMeDAS URL]_ |

Bundled `data/` folders contain the inputs used in the published analysis.

---

## 3. Usage

### 3.1 Requirements

- Python 3.9+
- `pip install -r requirements.txt`

### 3.2 Principles

1. Run each script from its region folder; inputs live in `./data/`.
2. Do not mix regions — workflows are independent.
3. Outputs (CSV, figures) are written to the working directory.
4. No random seeds; results depend on data and AIC-based model selection.

### 3.3 Run

**GBA mainland**

```bash
cd GBA
python contour_Estuary.py
```

**Hong Kong**

```bash
cd GBA/HongKong
python contour_HongKong.py
```

**San Francisco**

```bash
cd San_Francisco
python contour_SF.py
```

**Tokyo**

```bash
cd Tokyo
python contour_Tokyo.py
```

### 3.4 Customisation

| Parameter | Location | Purpose |
|---|---|---|
| `station_configs` | end of each script | Station pairs, plot limits |
| `enabled_distributions` | per station | Marginal candidates |
| `enabled_copulas` | per station (Estuary) | Copula candidates |
| `datum_offset` | SF configs | Vertical datum (m) |
| `offset_file` | Estuary / Tokyo | Datum correction series |
| `VERBOSE` | `contour_Estuary.py` | Verbose logging |

### 3.5 Outputs

Per tide station:

- `maximum_density_points_<STATION>.csv`
- `contour_axis_intersections_<STATION>.csv`
- Matplotlib figure (display or save manually)

### 3.6 Troubleshooting

| Issue | Check |
|---|---|
| `FileNotFoundError` | Input files present in `data/` |
| Empty compound events | Tide–rain overlap (see coverage figures) |
| GBA xlsx errors | Install `openpyxl` |
| SF missing rainfall | Rain gauge CSV files in `data/` |

---

## Repository layout

```
├── README.md
├── requirements.txt
├── docs/figures/
├── GBA/
│   ├── contour_Estuary.py
│   ├── data/
│   └── HongKong/
│       ├── contour_HongKong.py
│       └── data/
├── San_Francisco/
│   ├── contour_SF.py
│   └── data/
└── Tokyo/
    ├── contour_Tokyo.py
    └── data/
```

## License

Add your license file before public release (e.g. MIT).
