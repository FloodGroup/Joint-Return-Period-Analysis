"""Main program: run the five steps in order.

The program you run is a region's `run_analysis.py`:

    python gba/run_analysis.py --manifest MANIFEST --output OUTPUT
    python san_francisco/run_analysis.py --manifest MANIFEST --output OUTPUT
    python tokyo/run_analysis.py --manifest MANIFEST --output OUTPUT

Each of those scripts calls the five folders below. This module is the same
sequence in one place, used by `jrpa run` and `jrpa demo`.
"""

from jrpa.step01_prepare import step01_prepare
from jrpa.step02_events import step02_events
from jrpa.step03_annual import step03_annual
from jrpa.step04_copula import step04_copula
from jrpa.step05_export import step05_export


def run(manifest_path, output, settings=None, skip_compound=False):
    """Run steps 01 to 05 in order. Each region uses this same sequence."""
    analysis = step01_prepare(manifest_path, output, settings)
    step02_events(analysis)
    step03_annual(analysis)
    step04_copula(analysis, skip_compound=skip_compound)
    return step05_export(analysis)
