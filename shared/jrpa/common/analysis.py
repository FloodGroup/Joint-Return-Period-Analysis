"""In-memory products passed from step 01 through step 05."""

from dataclasses import dataclass, field
from pathlib import Path
from .config import Settings


@dataclass
class Analysis:
    """In-memory products passed through the five numbered processing steps."""

    manifest_path: Path
    manifest: dict
    settings: Settings
    out: Path
    inputs: list = field(default_factory=list)
    tide: dict = field(default_factory=dict)
    tide_quality: dict = field(default_factory=dict)
    water_maxima: dict = field(default_factory=dict)
    donors: list = field(default_factory=list)
    units: list = field(default_factory=list)
    sources: list = field(default_factory=list)
    annual_tide: dict = field(default_factory=dict)
    design_rows: list = field(default_factory=list)
    fitted: list = field(default_factory=list)
    completed: int = 1

    def require(self, step):
        if self.completed != step - 1:
            raise ValueError(f"Step {step:02d} requires completed step {step - 1:02d}")
