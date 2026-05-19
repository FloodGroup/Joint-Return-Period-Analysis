# Compound Flood Return-Period Analysis

Open-source scripts for estimating **bivariate rainfall–sea level return periods** in three megacity bay areas:

- **GBA_Pearl River Estuary** — Pearl River estuary (mainland China)
- **GBA_HongKong** — Hong Kong
- **San Francisco** — San Francisco Bay
- **Tokyo** — Tokyo Bay

---

## 1. Method

Each regional script follows the same steps:

1. **Pair stations** — Assign one or more rainfall gauges to each tide gauge (fixed pairs in GBA_Pearl River Estuary; priority-ordered lists in GBA_HongKong and Tokyo; San Francisco uses the six nearest gauges by distance).
2. **Merge and align** — Build a single daily rainfall series from paired gauges (first available date wins by priority); inner-join with daily maximum sea level on calendar date.
3. **Datum correction** — Convert daily maximum sea level to each region's unified vertical reference (see §1.1).
4. **Define compound events** — Flag days when 24-hour cumulative rainfall and daily maximum sea level both exceed their **95th-percentile** thresholds in the merged series.
5. **De-cluster** — Treat events within **3 days** of each other as one event.
6. **Fit models** — Select marginal distributions (AIC among Gamma, GEV, GP, logistic, normal, Gumbel, Nakagami) and copulas (AIC among Gaussian, Student-t, Clayton, Gumbel, Frank, Joe); estimate by MLE; assess marginals with the Kolmogorov–Smirnov test.
7. **Return periods** — For joint exceedance at levels (r, s):

$$
T(R>r,\ S>s)=\frac{1}{1-F_R(r)-F_S(s)+C(F_R(r),F_S(s))}
$$

where $F_R$, $F_S$ are marginal CDFs and $C$ is the fitted copula.
8. **Output** — Return-period contour plots; CSV tables of contour axis intersections and maximum-density points.

### 1.1 Vertical datum conversion

Tide-gauge sea levels are converted to a **regional unified vertical reference** before copula fitting. Rainfall is not transformed.

| Region | Unified reference |
|---|---|
| GBA_Pearl River Estuary | 1985 Chinese Height Datum (1985 中国国家高程基准) |
| GBA_HongKong | 1985 Chinese Height Datum (1985 中国国家高程基准) |
| San Francisco | North American Vertical Datum of 1988 (NAVD88) |
| Tokyo | Tokyo Peil (TP) |

**GBA_Pearl River Estuary** — two steps:

1. **Piecewise offset (Pearl River Datum)** — Original tide records are referenced to the Pearl River Datum (珠江基准面). Where a station’s vertical reference changed over time, a time-varying offset is applied to each daily maximum sea level. 

2. **Shift to 1985 datum** — After the piecewise step, all stations receive a uniform **+0.744 m** shift, equal to the height difference between the Pearl River Datum and the 1985 Chinese Height Datum.

**GBA_HongKong** — tide records are on Hong Kong Chart Datum (香港海图基准面, HKCD). All stations receive a uniform **−0.868 m** shift to the 1985 Chinese Height Datum (1985 datum is 0.868 m above HKCD).

**San Francisco Bay** — each tide gauge originally uses a different local datum. A station-specific constant offset converts all records to NAVD88. Negative offsets mean the station datum lies below NAVD88.

**Tokyo Bay** — same piecewise approach as GBA_Pearl River Estuary, with a station-specific offset series for each tide gauge. Corrected levels are on Tokyo Peil.

---

## 2. Data pairing and coverage

Coverage timelines are in [docs/figures/coverage_GBA.png](docs/figures/coverage_GBA.png), [docs/figures/coverage_San_Francisco.png](docs/figures/coverage_San_Francisco.png), and [docs/figures/coverage_Tokyo.png](docs/figures/coverage_Tokyo.png). Blue = tide observations; green = rainfall observations; red hatching = dates that rainfall actually supplies after priority merge and tide inner join. Percentage are the proportion of rainfall days at this station out of the total number of days selected for subsequent analysis.

### 2.1 GBA_Pearl River Estuary

| Tide gauge | Rainfall pairing | Tide observations | Analysis sample |
|---|---|---|---|
| Chiwan | Tiegang (100%) | 1965–2022 · 13,878 d | 1975–2022 · 10,226 d |
| Sanzao | Sanzao (100%) | 1965–2022 · 15,552 d | 1965–2022 · 15,552 d |
| Nei Lingding | Changjiang (100%) | 2010–2022 · 4,016 d | 2010–2022 · 4,016 d |

Pearl River Estuary rainfall workbooks treat blank cells as 0 mm daily rainfall. Tide levels are converted from Pearl River Datum to the 1985 Chinese Height Datum (§1.1).

### 2.2 GBA_HongKong

Tide levels are converted from Hong Kong Chart Datum (香港海图基准面) to the 1985 Chinese Height Datum (−0.868 m; §1.1).

| Tide gauge | Rainfall pairing | Tide observations | Analysis sample |
|---|---|---|---|
| Quarry Bay | Quarry Bay (27.5%) → Kai Tak (0.3%) → Shau Kei Wan (0.1%) → Hong Kong Observatory (72.1%) | 1960–2023 · 22,487 d | 1960–2023 · 22,487 d |
| Tsim Bei Tsui | Tsim Bei Tsui (94.6%) → Wetland Park (1.0%) → Lau Fau Shan (4.4%) | 1983–2023 · 12,847 d | 1985–2023 · 12,281 d |
| Shek Pik | Ngong Ping Fresh Water Reservoir (100%) | 1998–2023 · 9,073 d | 2006–2023 · 6,147 d |
| Tai Mo Wan | Sai Kung (Hong Kong Adventist College) (98.2%) → Shau Kei Wan (1.8%) | 1996–2023 · 9,076 d | 1996–2023 · 8,904 d |
| Tai Po Kau | Tai Po Wong Shiu Chi Secondary School (100%) | 1981–2023 · 15,248 d | 1985–2023 · 13,091 d |
| Waglan Island | Waglan Island (100%) | 1982–2023 · 8,980 d | 1989–2023 · 7,869 d |

### 2.3 San Francisco Bay

Six nearest rainfall gauges are assigned per tide station; gauges contributing less than 2% of the merged sample are omitted from analysis. Tide levels are converted to NAVD88 via station-specific offsets (§1.1).

| Tide gauge | Rainfall pairing | Tide observations | Analysis sample |
|---|---|---|---|
| San Francisco | San Francisco Downtown (83.0%) → Kentfield (10.3%) → Berkeley (6.6%) | 1898–2021 · 44,331 d | 1898–2021 · 43,259 d |
| Redwood City | Palo Alto (76.5%) → Newark (13.6%) → Woodside Fire Station (7.1%) → Fremont (2.7%) | 1974–2021 · 9,569 d | 1974–2021 · 9,569 d |
| Alameda | Oakland Museum (92.0%) → Oakland Metro Intl AP (5.3%) → San Francisco Downtown (2.6%) | 1976–2021 · 15,366 d | 1976–2021 · 15,366 d |
| Richmond | Richmond (90.8%) → Kentfield (9.2%) | 1996–2021 · 9,065 d | 1996–2021 · 9,065 d |
| Port Chicago | Concord Buchanan Field (53.5%) → Martinez Water Plant (46.3%) | 1979–2021 · 15,038 d | 1979–2021 · 15,038 d |

### 2.4 Tokyo Bay

Tide levels are converted to Tokyo Peil via station-specific piecewise offset series (§1.1).

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
| GBA_Pearl River Estuary | _[TODO: provider URL]_ | _[TODO: provider URL]_ |
| GBA_HongKong | _[TODO: HKO / CEDD URL]_ | _[TODO: HKO rainfall URL]_ |
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
2. Outputs (CSV, figures) are written to the working directory.
3. Results depend on data and AIC-based model selection.

### 3.3 Customisation

| Parameter | Location | Purpose |
|---|---|---|
| `station_configs` | end of each script | Station pairs, plot limits |
| `enabled_distributions` | per station | Marginal candidates |
| `enabled_copulas` | per station (Pearl River Estuary) | Copula candidates |
| `datum_offset` | SF configs | Station datum → NAVD88 shift (m) |
| `offset_file` | Pearl River Estuary / Tokyo | Piecewise datum correction series |
| `VERBOSE` | `contour_Estuary.py` | Verbose logging |

### 3.4 Outputs

Per tide station:

- `maximum_density_points_<STATION>.csv`
- `contour_axis_intersections_<STATION>.csv`
- Matplotlib figure (display or save manually)

---

## License

Add your license file before public release (e.g. MIT).
