"""Command-line entry point; private observations are always supplied explicitly."""

import argparse
import json
from pathlib import Path
from .common.config import demonstration_settings
from .common.io import write_json, write_csv


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="jrpa", description="Rainfall and coastal water-level frequency analysis"
    )
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="Analyze a private observation manifest")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--annual-only", action="store_true")
    demo = commands.add_parser("demo", help="Generate and analyze synthetic observations")
    demo.add_argument("--region", choices=["gba", "san_francisco", "tokyo"], default="tokyo")
    demo.add_argument("--output", type=Path, required=True)
    demo.add_argument("--years", type=int, default=16)
    normalize = commands.add_parser("normalize", help="Screen a private tabular provider export")
    normalize.add_argument("--input", type=Path, required=True)
    normalize.add_argument("--output", type=Path, required=True)
    normalize.add_argument("--variable", choices=["rain", "water"], required=True)
    normalize.add_argument(
        "--policy",
        choices=[
            "jma_rain",
            "mlit_rain",
            "sf_rain",
            "jma_tide",
            "gesla_tide",
            "jodc_tide",
            "verified",
        ],
        required=True,
    )
    normalize.add_argument("--time-column", default="time_utc")
    normalize.add_argument("--value-column")
    normalize.add_argument("--flag-column", default="flag")
    normalize.add_argument("--scale", type=float, default=1.0)
    normalize.add_argument("--daily", action="store_true")
    catalog = commands.add_parser(
        "catalog", help="Print public station identifiers and source priorities"
    )
    catalog.add_argument("--region", choices=["gba", "san_francisco", "tokyo"], required=True)
    cases = commands.add_parser(
        "cases", help="Generate scenario definitions without observation values"
    )
    cases.add_argument("--region", choices=["gba", "san_francisco", "tokyo"], required=True)
    cases.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "normalize":
        from .common.providers import normalize_export

        if args.input.resolve() == args.output.resolve():
            parser.error("Normalization output must not overwrite the provider export")
        frame, audit = normalize_export(
            args.input,
            args.policy,
            args.variable,
            args.time_column,
            args.value_column,
            args.flag_column,
            args.scale,
            args.daily,
        )
        write_csv(args.output, frame)
        write_json(args.output.with_suffix(".quality.json"), audit)
        print(json.dumps(audit, indent=2))
        return
    if args.command == "catalog":
        from .common.config import read_catalog

        print(json.dumps(read_catalog(args.region), indent=2))
        return
    if args.command == "cases":
        from .common.scenarios import case_definitions

        write_json(args.output, case_definitions(args.region))
        return
    from .main import run as analyze

    if args.command == "demo":
        from .common.synthetic import create_demo

        manifest = create_demo(args.output / "synthetic_inputs", args.region, args.years)
        result = analyze(manifest, args.output / "analysis", demonstration_settings(args.region))
    else:
        result = analyze(args.manifest, args.output, skip_compound=args.annual_only)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
