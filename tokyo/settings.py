"""Analysis parameters and station rules for Tokyo Bay."""

PARAMETERS = {
    "region": "tokyo",
    "datum": "TP",
    "timezone": "Asia/Tokyo",
    "event_clock": "hourly",
    "dry_gap_hours": 12,
    "decluster_hours": 72,
    "rainfall_quantile": 0.95,
    "minimum_water_coverage": 0.75,
    "minimum_events": 30,
    "minimum_effective_years": 10.0,
    "marginal_screen_reps": 49,
    "marginal_final_reps": 999,
    "joint_reps": 99,
    "annual_reps": 99,
    "seed": 202610052,
    "annual_seed": 20261015,
    "annual_coverage": 0.75,
    "minimum_annual_trial": 8,
    "minimum_annual_accepted": 15,
    "rainfall_factor_24h": 1.0,
    "kendall_samples": 1048576,
    "independent_calibration": True,
    "kendall_search": "global",
    "periods": (10, 50, 100, 200),
    "rain_peak_hour": 12.0,
    "synthetic": False,
}

TIDE_TRANSFERS = {"Kawasaki": "Tokyo-Harumi", "Kurihama": "Aburatsubo"}
