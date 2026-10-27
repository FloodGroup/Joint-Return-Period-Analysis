# Inputs and sources

Place observations and a `manifest.json` in a private directory outside the repository. Relative paths in the manifest resolve only inside that directory. Units are millimetres for rainfall and metres for water level.

## Canonical observation tables

| Series | Required columns | Clock and value definition |
|:--|:--|:--|
| Hourly rainfall | `time_utc`, `rain_mm` | UTC interval end; observed 1 h or intact 6 h accumulation |
| Hourly coastal water level | `time_utc`, `water_m` | UTC clock hour; target-datum observation |
| GBA rainfall | `date`, `rain_mm` | Fixed observation-day daily accumulation |
| GBA coastal water level | `date`, `water_m` | Local date; highest observed water level that day |

Missing values are empty fields, never zero. Timestamps must be unique and sorted by the reader. Do not disaggregate six-hour totals into invented hourly values. Provider quality flags, accumulation-period definitions, and datum corrections must be resolved explicitly during normalization.

## Manifest example

Example manifest for one Tokyo unit:

```json
{
  "region": "tokyo",
  "synthetic": false,
  "tides": {
    "Aburatsubo": {"path": "tides/aburatsubo.csv", "datum": "TP"}
  },
  "units": [{
    "unit": "M01",
    "gauge": "Aburatsubo",
    "rainfall_sources": [{
      "station_id": "JMA_A46136",
      "path": "rainfall/JMA_A46136.csv",
      "priority": 1,
      "resolution_h": 1
    }]
  }]
}
```

Supply every contributing source in the unit's declared order. `jrpa catalog --region tokyo` lists the study identifiers. The catalogs list contributing gauges; reports and maps omit unused sources.

## Epoch-specific vertical conversion

Instead of declaring an already normalized `datum`, a tide entry may provide `epochs`. Each epoch requires `start`, exclusive `end`, `offset_m`, `target_datum`, and an evidence `source`. An optional `epoch_timezone` defines naive epoch boundaries; Tokyo uses JST. The transformation is

$$H_{\mathrm{target}}(t)=H_{\mathrm{source}}(t)+\Delta z_e,\quad t\in[t_{e,0},t_{e,1}).$$

Unknown offsets remain missing. Overlapping epochs or mismatched target datums cause an error. A normalized tide series must not receive the conversion a second time.

## Quality-control fields

Rainfall and tide entries accept `exclusions`, each containing an inclusive `start`, `end`, and a reason. These affect the common normalized source. `annual_exclusions` apply only to the full-record annual calculation, preserving the event population used by the compound analysis. Values removed by a rule remain auditable in the private output provenance.

For San Francisco, the annual quality controls include the documented Alviso terminal ramp, reviewed Mare Island/Martinez intervals, and isolated rainfall spikes. Provider-specific flag meanings must be interpreted before canonical input is created.

## GBA-specific inputs

Mainland S1 units additionally require a `rain_location` pair `[latitude, longitude]` and a `regional_rainfall` list. Each regional station entry contains `station_id`, `path`, and `location`. The code computes annual maxima and uses only long-record stations within 60 km of the representative location, subject to the declared heterogeneity screening.

Set each GBA unit's `event_span` to its declared `[start_date, end_date]` analysis clock. Missing tide dates remain on that clock: positive rainfall still contributes to the R95 threshold population, while paired exposure counts only dates with both variables observed. The full rainfall record remains available separately for S1.

The tide manifest includes Hong Kong reference gauges for storm-day checks and both short stations' transfer donors. `supplemental_peaks` points to a private JSON table of independently sourced peak observations. Each row contains `station_id`, `year`, `water_m`, and an evidence `source`. Keep this observation table outside the repository.

## Normalize and retain provenance

Tabular provider exports can be screened with the normalization command:

```sh
jrpa normalize --input ../private/raw/jma-export.csv --output ../private/tokyo/rainfall/JMA_A46136.csv --variable rain --policy jma_rain
```

The default source columns are `time_utc`, `precip_mm` or `water_level_m`, and `flag`. Explicit column options and `--scale` handle tabular export variants and units. `--daily` uses date labels. The output retains rejected values as missing and stores a quality audit. It does not apply vertical offsets.

JMA rain keeps quality codes beginning with 8 or 5. San Francisco rain rejects negative values and A flags. MLIT rain keeps confirmed/provisional export flags. JMA tide retains the official reanalysis R code; JODC tide rejects its anomalous R flag; GESLA tide keeps good code 1. These meanings must not be interchanged. The `verified` policy is only for an independently screened export.

1. Acquire provider records and station/datum documentation independently.
2. Preserve raw exports in the private archive.
3. Convert units, validate quality codes, resolve duplicate timestamps, and document accumulation intervals.
4. Apply documented datum epochs or declare the verified target datum.
5. Run `jrpa run --manifest PATH --output PATH`.

Provider-native exports must be converted to the canonical columns before analysis. Preserve missing values and documented datum offsets. Input hashes and source contributions are recorded with each run.

## Observation providers

| Region | Variable | Provider and record family |
|:--|:--|:--|
| Greater Bay Area | Mainland coastal water level | Pearl River Basin Hydrological Yearbook, daily highest coastal water level |
| Greater Bay Area | Mainland rainfall | [China Meteorological Administration](https://data.cma.cn/), hydrological-yearbook stations, and [NOAA NCEI](https://www.ncei.noaa.gov/) long-record station products |
| Greater Bay Area | Hong Kong rain and tide | [Hong Kong Observatory](https://www.hko.gov.hk/en/cis/climat.htm), daily rainfall and daily highest high water |
| San Francisco Bay | Coastal water level | [NOAA CO-OPS](https://tidesandcurrents.noaa.gov/), [GESLA](https://gesla787883612.wordpress.com/), and [USGS Water Data](https://waterdata.usgs.gov/) |
| San Francisco Bay | Rainfall | [California Data Exchange Center](https://cdec.water.ca.gov/) and [NOAA NCEI](https://www.ncei.noaa.gov/) station observations |
| Tokyo Bay | Coastal water level | [Japan Oceanographic Data Center](https://www.jodc.go.jp/), [JMA](https://www.data.jma.go.jp/kaiyou/db/tide/sea_lev_var/), and GESLA station archives |
| Tokyo Bay | Rainfall | [JMA historical weather observations](https://www.data.jma.go.jp/stats/etrn/index.php), including AMeDAS station identifiers |

Station identifiers preserve provider namespaces. For example, `NOAA_9414523`, `USGS_372750122012701`, `CDEC_CMG`, `JMA_A46136`, and `GESLA_JODC-GIAJ_GS01` identify sources without redistributing their observations. Provider aliases such as GBA `SZ`, `DJZ`, and `NLT` are defined in the regional documentation.


GBA regional rainfall reference stations include Sanzao, Changjiang, Zhongda, Henggang, Tiegang, Nansha, Guangzhou (`GHCND_CHM00059287`) and Baoan (`ISD_594930`). Supply verified coordinates and annual samples in the private regional station table.

Record each product's download date, version, station identifier, time zone, units and quality-code interpretation. Station identifiers alone do not fix a provider snapshot. The software license does not grant rights to provider observations.

## Reviewed annual samples

An optional `annual_maxima` path on a unit or tide entry supplies reviewed annual observations. Its columns are `year`, `maximum` and Boolean `used`. Years must be unique and accepted maxima finite. This input preserves documented spike exclusions, partial-year decisions and boundary-window selection. Fitted design levels are not accepted as annual observations.

## Study-area figures

The regional pages include analysis polygons, station positions, rainfall priorities and observation availability. Each region's `stations.json` contains station names and coordinates. The drawn figures are in that region's `figures/` folder. Source GIS and the drawing program are not part of this release.

Roles are `rain`, `tide` and `paired`. Fractions use the entire calendar month as denominator. Rainfall rows represent selected observations from contributing gauges. A `station_id` of `merged_rain` supplies the unit's combined rainfall row.

First-use and last-use dates bound actual contributions to the archived compound-input record. Gaps are shown in the coverage figures. These dates do not describe the complete provider archive or the independently screened annual-maximum sample. GBA gauges without a documented name retain their provider identifier.
