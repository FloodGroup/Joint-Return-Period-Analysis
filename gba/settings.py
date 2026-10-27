"""Analysis parameters and station rules for Greater Bay Area."""

PARAMETERS = {
    "region": "gba",
    "datum": "CN1985",
    "timezone": "Asia/Hong_Kong",
    "event_clock": "daily",
    "dry_gap_hours": 12,
    "decluster_hours": 72,
    "rainfall_quantile": 0.95,
    "minimum_water_coverage": 0.75,
    "minimum_events": 30,
    "minimum_effective_years": 10.0,
    "marginal_screen_reps": 49,
    "marginal_final_reps": 49,
    "joint_reps": 49,
    "annual_reps": 99,
    "seed": 20261011,
    "annual_seed": 20261015,
    "annual_coverage": 0.75,
    "minimum_annual_trial": 8,
    "minimum_annual_accepted": 15,
    "rainfall_factor_24h": 1.13,
    "kendall_samples": 262144,
    "independent_calibration": False,
    "kendall_search": "daily_grid",
    "periods": (10, 50, 100, 200),
    "rain_peak_hour": 12.0,
    "synthetic": False,
}

TIDE_TRANSFERS = {"DJZ": "SZ", "NLT": "CW"}

GBA_EVENT_EXCLUSIONS = {"01": ["2014-09-16"], "02": ["2014-09-16"], "05": ["2017-08-23"]}

# Regional index-rainfall analysis.
REGIONAL_RADIUS_KM = 60
REGIONAL_H1_LIMIT = 2
REGIONAL_REPLICATES = 500
REGIONAL_SEED = 7

# Analysis datum to model datum, in metres.
GBA_MAINLAND_MODEL_OFFSET = -0.744
GBA_HK_MODEL_OFFSET = -0.573
