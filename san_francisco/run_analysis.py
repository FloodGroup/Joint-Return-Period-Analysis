"""Run san_francisco rainfall and coastal water-level frequency analysis."""

import argparse
import json
from pathlib import Path
from jrpa.common.config import Settings
from jrpa.san_francisco.settings import PARAMETERS
from jrpa.step01_prepare import step01_prepare
from jrpa.step02_events import step02_events
from jrpa.step03_annual import step03_annual
from jrpa.step04_copula import step04_copula
from jrpa.step05_export import step05_export


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--annual-only", action="store_true")
    args = parser.parse_args()
    settings = Settings(**PARAMETERS)

    # 01. Read observations; apply quality rules and vertical references.
    analysis = step01_prepare(args.manifest, args.output, settings)
    # 02. Select rainfall sources and construct independent event pairs.
    step02_events(analysis)
    # 03. Estimate S1 rainfall and S2 water-level frequency curves.
    step03_annual(analysis)
    # 04. Fit event margins and copula dependence.
    step04_copula(analysis, skip_compound=args.annual_only)
    # 05. Calculate compound design values and export the report.
    summary = step05_export(analysis)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
