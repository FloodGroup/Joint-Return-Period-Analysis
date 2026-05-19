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
3. **Datum correction** — Convert daily maximum sea level between vertical datums to the regional target reference (see §1.1).
4. **Define extreme events** — Flag days when 24-hour cumulative rainfall and daily maximum sea level both exceed their **95th-percentile** thresholds in the merged series.
5. **De-cluster** — Treat events within **3 days** of each other as one event.
6. **Fit models** — Select marginal distributions (AIC among Gamma, GEV, GP, logistic, normal, Gumbel, Nakagami) and copulas (AIC among Gaussian, Student-t, Clayton, Gumbel, Frank, Joe); estimate by MLE; assess marginals with the Kolmogorov–Smirnov test.
7. **Return periods** — For joint exceedance at levels (r, s):

$$
T(R>r,\ S>s)=\frac{1}{1-F_R(r)-F_S(s)+C(F_R(r),F_S(s))}
$$

where $F_R$, $F_S$ are marginal CDFs and $C$ is the fitted copula.
8. **Output** — Return-period contour plots; CSV tables of contour axis intersections and maximum-density points.

### 1.1 Vertical datum conversion

Tide-gauge records were originally collected on different vertical datums. Before copula fitting, sea levels are converted from their source datum to a common target datum in each region.

| Region | Source datum | Target datum |
|---|---|---|
| GBA_Pearl River Estuary | Pearl River Datum (珠江基准面), with time-varying adjustments | 1985 Chinese Height Datum (1985 中国国家高程基准) |
| GBA_HongKong | Hong Kong Chart Datum (香港海图基准面) | 1985 Chinese Height Datum (1985 中国国家高程基准) |
| San Francisco | Station-specific local datums | North American Vertical Datum of 1988 (NAVD88) |
| Tokyo | Station-specific local datums, with time-varying adjustments | Tokyo Peil (TP) |

**GBA_Pearl River Estuary** — two steps:

1. **Piecewise offset (Pearl River Datum)** — Original tide records are referenced to the Pearl River Datum (珠江基准面). Where a station’s vertical reference changed over time, a time-varying offset is applied to each daily maximum sea level.

2. **Shift to 1985 datum** — After the piecewise step, all stations receive a uniform **+0.744 m** shift, equal to the height difference between the Pearl River Datum and the 1985 Chinese Height Datum.

**GBA_HongKong** — tide records are on Hong Kong Chart Datum (香港海图基准面, HKCD). All stations receive a uniform **−0.868 m** shift to the 1985 Chinese Height Datum (1985 datum is 0.868 m above HKCD).

**San Francisco Bay** — each tide gauge originally uses a different local datum. A station-specific constant offset converts all records to NAVD88. Negative offsets mean the station datum lies below NAVD88.

**Tokyo Bay** — same piecewise approach as GBA_Pearl River Estuary, with a station-specific offset series for each tide gauge. Corrected levels are on Tokyo Peil.

---

## 2. Data pairing and coverage

Each figure covers **one region**. Within a figure, **each subplot corresponds to one tide gauge used in the analysis** and shows (i) the tide record span, (ii) the observation spans of rainfall gauges paired with that tide gauge, and (iii) the matched period used after priority merging and tide–rainfall inner join. **Blue** bars mark contiguous tide-observation periods; **green** bars mark rainfall records paired with that gauge; **red hatching** marks the dates that a rainfall gauge actually contributes to the merged series. The percentage after each rainfall station name is its share of the analysable sample at that tide gauge.

### 2.1 GBA

![GBA tide–rainfall coverage](docs/figures/coverage_GBA.png)

**GBA_Pearl River Estuary** (Chiwan, Sanzao, Nei Lingding) uses fixed one-to-one tide–rainfall pairs. Chiwan is limited by rainfall from 1975 onward (10,226 analysis days); Sanzao spans 1965–2022 (15,552 days); Nei Lingding covers 2010–2022 (4,016 days). **GBA_HongKong** (six tide gauges) fills gaps with priority-ordered rainfall lists—for example, Quarry Bay draws mainly from Hong Kong Observatory (72%) and Quarry Bay (28%). Analysis periods range from about 6,000 days (Shek Pik) to 22,500 days (Quarry Bay).

### 2.2 San Francisco Bay

![San Francisco tide–rainfall coverage](docs/figures/coverage_San_Francisco.png)

Five tide gauges; each is matched to the six nearest rainfall gauges, keeping sources that contribute at least 2% of the merged sample. San Francisco has the longest record (1898–2021, 43,259 analysis days). The other four stations span 1974–2021 with 9,000–15,400 analysis days each.

### 2.3 Tokyo Bay

![Tokyo Bay tide–rainfall coverage](docs/figures/coverage_Tokyo.png)

Six tide gauges with manually assigned rainfall priority lists. Most stations yield 15,000–20,000 analysis days from the 1960s–2010s; Yokosuka and Aburatsubo are shorter from 1976 onward because rainfall records start later. Kawasaki ends in 1996 (8,771 days).

### 2.4 Data sources

| Region | Tide data | Rainfall data |
|---|---|---|
| GBA_Pearl River Estuary | Pearl River Basin Hydrological Yearbook (珠江流域水文年鉴) | Pearl River Basin Hydrological Yearbook (珠江流域水文年鉴) |
| GBA_HongKong | [Hong Kong Observatory — Climatological Information Services](https://my.weather.gov.hk/en/cis/climat.htm) | [Hong Kong Observatory — Climatological Information Services](https://my.weather.gov.hk/en/cis/climat.htm) |
| San Francisco | [GESLA-4](https://gesla787883612.wordpress.com/downloads/) | [NOAA NCEI — Daily Observational Data](https://www.ncei.noaa.gov/maps/daily/) |
| Tokyo | [GESLA-4](https://gesla787883612.wordpress.com/downloads/) | [JMA — Past Weather Data Search](https://www.data.jma.go.jp/stats/etrn/index.php) |

Bundled `data/` folders contain the station files used in the published analysis.

---

## 3. Usage

### 3.1 Requirements

- Python 3.9+
- `pip install -r requirements.txt`

### 3.2 Principles

1. **Run the script for the target bay area only** — Do not mix regions: input formats and vertical-datum conversion all differ (§1.1).
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
