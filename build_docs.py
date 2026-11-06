"""Build HTML documentation inside each study-area folder."""

import argparse
from functools import cache
import hashlib
import io
from pathlib import Path
import re
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
:root{--md-primary-fg-color:#244957;--md-accent-fg-color:#087f86;--md-default-fg-color:#24343d;--md-default-fg-color--light:#52636d}
body{font-family:"Segoe UI",Arial,sans-serif}
.md-grid{max-width:80rem}.md-sidebar{width:10.8rem}.md-sidebar.md-sidebar--secondary{display:none!important}
.md-typeset{font-size:.85rem;line-height:1.75}
.md-typeset h1{font-size:1.9rem;line-height:1.2;font-weight:650;letter-spacing:-.025em;color:#203844;margin:0 0 1.2rem}
.md-typeset h2{font-size:1.3rem;font-weight:650;line-height:1.35;letter-spacing:-.015em;margin:2.4rem 0 1rem;padding-top:.7rem;border-top:1px solid #dae3e7;scroll-margin-top:4rem}
.md-typeset h3{font-size:1.05rem;font-weight:650;line-height:1.45;margin:1.8rem 0 .8rem;scroll-margin-top:4rem}
.md-typeset h3[id^="unit-"]{font-size:1.3rem;border-top:2px solid #9fb9c1;padding-top:1.3rem;margin-top:3rem}
.md-typeset h4{font-size:.85rem;font-weight:600;color:#425b66;margin:1.3rem 0 .6rem}
.md-typeset p{margin:.7rem 0 1rem}.md-typeset li{margin-bottom:.45rem}
.md-typeset .md-typeset__table{display:block;min-width:0}
.md-typeset table:not([class]){display:table;width:100%;font-size:.77rem;line-height:1.6;border:1px solid #dbe4e8;border-radius:4px;font-variant-numeric:tabular-nums}
.md-typeset table:not([class]) th{background:#f0f5f6;font-weight:650;min-width:0;padding:.7rem .8rem;border-bottom:1px solid #cbd9df;vertical-align:bottom}
.md-typeset table:not([class]) td{padding:.6rem .8rem;vertical-align:top;border-top:1px solid #e4eaed}
.md-typeset table:not([class]) tbody tr:nth-child(even){background:#fafcfc}
.md-typeset table:not([class]) td code{white-space:nowrap;font-size:.72rem}
.md-typeset .md-typeset__scrollwrap{margin:1rem 0 1.5rem}
.md-typeset .arithmatex{overflow-x:auto;padding:.8rem .5rem;margin:1.2rem 0;font-size:1em}
.md-typeset .katex-display{margin:.4em 0}.md-typeset code{font-size:.8em}
.md-typeset pre{margin:1rem 0 1.5rem;border:1px solid #dbe4e8;border-radius:4px}
.md-typeset pre>code{font-size:.75rem;line-height:1.7;padding:1rem}
.md-typeset img{display:block;max-width:100%;height:auto;margin:0 auto}
.md-typeset p:has(>a>img){margin:1.1rem 0 .4rem;padding:.4rem;background:white;border:1px solid #e2e9ec;border-radius:4px}
.md-typeset p:has(>a>img)+p:has(>em:only-child){font-size:.75rem;line-height:1.6;color:#52636d;margin:.4rem 0 1.6rem}
.md-typeset p:has(>a>img)+p>em{font-style:normal}
.md-header{box-shadow:none;border-bottom:1px solid #d5e0e5}
.md-header__topic:first-child{display:none}.md-header__topic:last-child{position:static;opacity:1;transform:none;transition:none}
.md-nav{font-size:.72rem;line-height:1.5}.md-nav__link{margin-top:.5rem}
.md-typeset strong{font-weight:650}.md-typeset a:has(>img){cursor:zoom-in;display:block}
@media(min-width:76.25em){.md-sidebar--secondary{display:none}.md-content{max-width:61rem;margin:0 auto;width:calc(100% - 11rem)}.md-content__inner{margin:0 1.8rem 3rem}}
@media(max-width:76.24em){.md-sidebar--secondary{display:none}.md-content{min-width:0}.md-content__inner{margin:0 1rem 2rem}}
@media(max-width:44.99em){.md-typeset{font-size:.82rem}.md-typeset h1{font-size:1.55rem}.md-typeset h2{font-size:1.2rem}.md-typeset table:not([class]){font-size:.73rem}.md-typeset table:not([class]) th,.md-typeset table:not([class]) td{padding:.55rem .65rem}.md-typeset p:has(>a>img){padding:0}.md-content__inner{margin:0 .8rem 2rem}}
@media print{.md-sidebar,.md-header,.md-tabs,.md-footer{display:none}.md-content{max-width:none!important;width:100%!important}.md-typeset h2,.md-typeset h3{break-after:avoid}.md-typeset img,.md-typeset table{break-inside:avoid}.md-typeset p:has(>a>img){border:0}.md-typeset .md-typeset__scrollwrap{overflow:visible}}
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


def documentation_text(root, source, region):
    """Resolve repository links for the five pages in each offline site."""
    folder = root / REGION_FOLDERS[region]
    pages = {
        folder / "README.md": "index.md",
        root / "README.md": "overview.md",
        root / "shared/README.md": "workflow.md",
        root / "shared/methods.md": "methods.md",
        root / "shared/inputs.md": "inputs.md",
    }

    def resolve(match):
        href = match.group(1)
        if href.startswith(("#", "https://", "http://", "mailto:")):
            return match.group(0)
        path, separator, fragment = href.partition("#")
        target = (source.parent / path).resolve()
        suffix = separator + fragment if separator else ""
        if target in pages:
            link = pages[target]
        elif target.is_relative_to(folder):
            link = target.relative_to(folder).as_posix()
        else:
            link = (
                "https://github.com/FloodGroup/Joint-Return-Period-Analysis/blob/main/"
                + target.relative_to(root).as_posix()
            )
        return "](" + link + suffix + ")"

    return re.sub(r"\]\(([^)]+)\)", resolve, source.read_text(encoding="utf8"))


def build_region(root, region):
    folder = root / REGION_FOLDERS[region]
    with tempfile.TemporaryDirectory(prefix="jrpa-docs-") as temporary:
        docs = Path(temporary) / "docs"
        docs.mkdir()
        for source, destination in (
            (folder / "README.md", "index.md"),
            (root / "README.md", "overview.md"),
            (root / "shared/README.md", "workflow.md"),
            (root / "shared/methods.md", "methods.md"),
            (root / "shared/inputs.md", "inputs.md"),
        ):
            (docs / destination).write_text(
                documentation_text(root, source, region), encoding="utf8"
            )
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
                {"Analysis workflow": "workflow.md"},
                {"Methods": "methods.md"},
                {"Inputs and sources": "inputs.md"},
                {"Repository overview": "overview.md"},
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
