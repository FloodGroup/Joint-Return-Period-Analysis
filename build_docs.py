"""Build HTML documentation inside each study-area folder."""

import argparse
from pathlib import Path
import shutil
import tempfile
from mkdocs.commands.build import build
from mkdocs.config import load_config

REGIONS = {
    "gba": "Greater Bay Area",
    "san_francisco": "San Francisco Bay",
    "tokyo": "Tokyo Bay",
}
REGION_FOLDERS = {"gba": "GBA", "san_francisco": "SanF", "tokyo": "Tokyo"}


def build_region(root, region):
    folder = root / REGION_FOLDERS[region]
    with tempfile.TemporaryDirectory(prefix="jrpa-docs-") as temporary:
        docs = Path(temporary) / "docs"
        docs.mkdir()
        text = (folder / "README.md").read_text(encoding="utf8")
        text = text.replace("../shared/methods.md", "methods.md")
        text = text.replace("../shared/inputs.md", "inputs.md")
        (docs / "index.md").write_text(text, encoding="utf8")
        for name in ("methods.md", "inputs.md"):
            shutil.copy2(root / "shared" / name, docs / name)
        for name in ("settings.py", "stations.json", "run_analysis.py"):
            shutil.copy2(folder / name, docs / name)
        shutil.copytree(root / "shared/assets", docs / "assets")
        shutil.copytree(folder / "figures", docs / "figures")
        config = load_config(
            config_file=str(root / "mkdocs.yml"),
            docs_dir=str(docs),
            site_dir=str(folder / "site"),
            site_name=REGIONS[region] + " | Joint Return-Period Analysis",
            nav=[
                {REGIONS[region]: "index.md"},
                {"Methods": "methods.md"},
                {"Inputs and sources": "inputs.md"},
            ],
            strict=True,
        )
        build(config)
    print(folder / "site/index.html")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", choices=REGIONS)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    for region in [args.region] if args.region else REGIONS:
        build_region(root, region)


if __name__ == "__main__":
    main()
