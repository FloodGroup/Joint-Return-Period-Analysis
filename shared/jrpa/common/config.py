"""Explicit, serializable settings for regional inference."""

from dataclasses import asdict, dataclass, replace
from importlib import import_module
from importlib.resources import files
import json


@dataclass(frozen=True)
class Settings:
    region: str
    datum: str
    timezone: str
    event_clock: str = "hourly"
    dry_gap_hours: int = 12
    decluster_hours: int = 72
    rainfall_quantile: float = 0.95
    minimum_water_coverage: float = 0.75
    minimum_events: int = 30
    minimum_effective_years: float = 10.0
    marginal_screen_reps: int = 49
    marginal_final_reps: int = 999
    joint_reps: int = 99
    annual_reps: int = 99
    seed: int = 202610052
    annual_seed: int = 20261015
    annual_coverage: float = 0.75
    minimum_annual_trial: int = 8
    minimum_annual_accepted: int = 15
    rainfall_factor_24h: float = 1.0
    kendall_samples: int = 1048576
    independent_calibration: bool = True
    kendall_search: str = "global"
    periods: tuple[int, ...] = (10, 50, 100, 200)
    rain_peak_hour: float = 12.0
    synthetic: bool = False

    def __post_init__(self):
        if (
            self.dry_gap_hours != 12
            or self.minimum_water_coverage != 0.75
            or self.decluster_hours != 72
        ):
            raise ValueError(
                "The study event kernel requires 12 h dry gaps, 75% water coverage, and 72 h declustering"
            )
        if self.rainfall_quantile != 0.95:
            raise ValueError("The study threshold is fixed at the 95th percentile")

    def as_dict(self):
        return asdict(self)


def region_module(region):
    """Load the settings stored alongside each region's entry scripts."""
    if region not in {"gba", "san_francisco", "tokyo"}:
        raise ValueError(f"Unknown region: {region}")
    return import_module(f"jrpa.{region}.settings")


def settings_for(region):
    return Settings(**region_module(region).PARAMETERS)


def tide_transfers(region):
    return region_module(region).TIDE_TRANSFERS


def read_catalog(region):
    region_module(region)
    return json.loads(files(f"jrpa.{region}").joinpath("stations.json").read_text(encoding="utf8"))


def demonstration_settings(region):
    """Small deterministic workload, explicitly ineligible as a study result."""
    return replace(
        settings_for(region),
        marginal_screen_reps=9,
        marginal_final_reps=19,
        joint_reps=9,
        annual_reps=19,
        kendall_samples=65536,
        independent_calibration=False,
        periods=(10, 50),
        synthetic=True,
    )
