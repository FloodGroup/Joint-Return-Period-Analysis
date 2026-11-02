"""Build HTML documentation inside each study-area folder."""

import argparse
from functools import cache
import hashlib
import io
from pathlib import Path
import shutil
import tempfile
import tarfile
from urllib.request import urlopen
from mkdocs.commands.build import build
from mkdocs.config import load_config

REGIONS = {
    "gba": "Greater Bay Area",
    "san_francisco": "San Francisco Bay",
    "tokyo": "Tokyo Bay",
}
REGION_FOLDERS = {"gba": "GBA", "san_francisco": "SanF", "tokyo": "Tokyo"}


KATEX_VERSION = "0.16.22"
KATEX_SHA256 = "e9e0d167db3175481cbadaff38e8d90b130f6a3ddb451a47e43c577fd511f365"
STYLES = """
:root{--md-primary-fg-color:#244957;--md-accent-fg-color:#168c98;--md-default-fg-color:#203844}
body{font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.md-typeset h1{font-weight:650;letter-spacing:-.025em;color:#203844}
.md-typeset h2{font-weight:600;letter-spacing:-.015em}
.md-typeset table:not([class]){font-size:.72rem;border-color:#e5eaed}
.md-typeset table:not([class]) th{background:#f2f6f7;font-weight:600;min-width:5rem}
.md-typeset .arithmatex{overflow-x:auto;padding:.4rem 0}
.md-typeset img{max-width:100%;height:auto}
.md-header{box-shadow:none;border-bottom:1px solid #d5e0e5}
.md-typeset strong{font-weight:650}

.md-typeset a:has(>img){cursor:zoom-in}
"""
MATH = r"""
document.addEventListener("DOMContentLoaded",function(){
  if(typeof renderMathInElement!=="function")return;
  renderMathInElement(document.body,{
    delimiters:[{left:"$$",right:"$$",display:true},{left:"\\[",right:"\\]",display:true},
                {left:"\\(",right:"\\)",display:false},{left:"$",right:"$",display:false}],
    throwOnError:true,strict:"error",trust:false
  });
});
"""


@cache
def katex_archive():
    """Fetch the pinned equation renderer once per documentation build."""
    url = f"https://registry.npmjs.org/katex/-/katex-{KATEX_VERSION}.tgz"
    with urlopen(url, timeout=60) as response:
        data = response.read()
    if hashlib.sha256(data).hexdigest() != KATEX_SHA256:
        raise ValueError("KaTeX archive checksum mismatch")
    return data


def create_assets(destination):
    """Generate documentation styles and renderer in the temporary build tree."""
    destination.mkdir()
    (destination / "styles.css").write_text(STYLES, encoding="utf8")
    (destination / "math.js").write_text(MATH, encoding="utf8")
    vendor = destination / "vendor/katex"
    with tarfile.open(fileobj=io.BytesIO(katex_archive()), mode="r:gz") as archive:
        for item in archive.getmembers():
            prefix = "package/dist/"
            if not item.isfile() or not item.name.startswith(prefix):
                continue
            relative = Path(item.name[len(prefix) :])
            if relative.name not in {
                "katex.min.js",
                "katex.min.css",
                "auto-render.min.js",
            } and not (
                relative.parts[0] == "fonts" and relative.suffix in {".woff2", ".woff", ".ttf"}
            ):
                continue
            # The upstream auto-render script lives under dist/contrib/.
            target = vendor / (relative.name if relative.name == "auto-render.min.js" else relative)
            if not target.resolve().is_relative_to(vendor.resolve()):
                raise ValueError("Invalid path in KaTeX archive")
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.extractfile(item) as source:
                target.write_bytes(source.read())
        license_file = archive.extractfile("package/LICENSE")
        with license_file:
            (vendor / "LICENSE").write_bytes(license_file.read())


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
        create_assets(docs / "assets")
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
