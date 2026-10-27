"""Offline English reports with consistent, vector-first scientific figures."""

from html import escape
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from jrpa.common.plotting import COLORS, PERIOD_COLORS, save, style
from jrpa.step04_copula.marginals import Margin

CSS = """
:root{--ink:#203844;--muted:#71818a;--teal:#168c98;--line:#e5eaed}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;color:var(--ink);font:15px/1.65 system-ui,sans-serif;background:#f5f7f8}
main{max-width:1200px;margin:40px auto;padding:0 28px}header{padding:32px 0;border-top:4px solid var(--teal)}
h1{font-size:clamp(28px,4vw,44px);line-height:1.15;letter-spacing:-.025em;margin:12px 0}h2{font-size:23px;margin:0 0 14px}
.eyebrow{font-size:11px;font-weight:700;letter-spacing:.15em;color:var(--teal);text-transform:uppercase}.muted,small{color:var(--muted)}
section{background:white;border:1px solid var(--line);padding:26px;margin:20px 0;border-radius:7px}.notice{background:#fff5e9;border-left:3px solid #ba7143;padding:14px 18px}
nav{display:flex;gap:10px;flex-wrap:wrap}a{color:#277383;text-decoration:none}a:hover{text-decoration:underline}
nav a{padding:6px 12px;border:1px solid var(--line);background:white;border-radius:4px}figure{margin:22px 0}img{width:100%;height:auto}figcaption{font-size:12px;color:var(--muted)}
table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}th,td{padding:9px 12px;border-bottom:1px solid var(--line);text-align:right}th:first-child,td:first-child{text-align:left}th{color:var(--muted);font-weight:600}thead{background:#f7f9fa}.scroll{overflow:auto}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px}.stat{background:white;border:1px solid var(--line);padding:15px 20px;border-radius:5px}.stat b{display:block;font-size:25px}
code{font-size:.9em;background:#edf3f5;padding:2px 4px;border-radius:3px}footer{padding:30px 0;color:var(--muted);font-size:12px}
"""


def unit_figures(root, meta, synthetic=False):
    style()
    root = Path(root)
    data = pd.read_csv(root / "pot_events.csv")
    curves = pd.read_csv(root / "kendall_curves.csv")
    points = pd.read_csv(root / "kendall_design.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.5), layout="constrained")
    for j, (variable, col, label) in enumerate(
        [
            ("R", "rain_mm", "Event rainfall (mm)"),
            ("H", "water_max_m", f"Coastal water level (m, {meta['datum']})"),
        ]
    ):
        x = np.sort(data[col].to_numpy(float))
        f = (np.arange(len(x)) + 0.5) / len(x)
        ax = axes[j]
        margin = Margin(meta["margins"][variable])
        ax.scatter(
            x,
            f,
            s=13,
            color=COLORS["sample"],
            label="Generated event sample" if synthetic else "Observed event sample",
        )
        grid = np.linspace(float(x.min()), float(max(x.max(), margin.ppf(0.995))), 240)
        ax.plot(
            grid,
            margin.cdf(grid),
            color=COLORS["rain" if variable == "R" else "water"],
            label="Fitted margin",
        )
        ax.set(xlabel=label, ylabel="Cumulative probability", ylim=(0, 1.025))
        ax.set_title(f"{'ab'[j]}  Marginal distribution", loc="left")
        ax.grid(alpha=0.8)
        ax.legend(fontsize=7, loc="lower right")
    ax = axes[2]
    factor = meta["settings"]["rainfall_factor_24h"]
    ax.scatter(
        data.rain_mm * factor, data.water_max_m, s=12, color=COLORS["sample"], alpha=0.75, zorder=1
    )
    for color, (period, part) in zip(PERIOD_COLORS, curves.groupby("return_period_years")):
        ax.plot(part.R_mm, part.H_m, color=color, label=f"{int(period)} yr")
        p = points.loc[points.return_period_years.eq(period)].iloc[0]
        ax.scatter(
            p.R_mm,
            p.H_m,
            color=color,
            marker="*",
            s=110,
            edgecolors="white",
            linewidths=0.5,
            zorder=5,
        )
    ax.set(xlabel="24 h rainfall (mm)", ylabel=f"Coastal water level (m, {meta['datum']})")
    ax.set_title("c  Kendall design", loc="left")
    ax.grid(alpha=0.8)
    ax.legend(fontsize=7)
    fig.suptitle(
        f"{meta['unit']} | "
        + ("Synthetic demonstration" if synthetic else "Rainfall-coastal water-level analysis"),
        x=0.01,
        ha="left",
        fontsize=13,
    )
    save(fig, root / "figures/model_diagnostics")


def render_report(output, settings, design_rows, models, sources):
    out = Path(output)
    title = {"gba": "Greater Bay Area", "san_francisco": "San Francisco Bay", "tokyo": "Tokyo Bay"}[
        settings.region
    ]
    frame = pd.DataFrame(design_rows)
    sections = []
    for meta in models:
        folder = out / "units" / meta["unit"]
        unit_figures(folder, meta, settings.synthetic)
        annual_figure(out, meta, settings)
        data = frame.loc[frame.unit.eq(meta["unit"])].copy()
        s2 = frame.loc[frame.scenario.eq("S2") & frame.gauge.eq(meta["gauge"])].copy()
        data = pd.concat([data, s2], ignore_index=True)
        columns = ["scenario", "return_period_years", "R_mm", "H_m", "H_model_m", "diagnostic_only"]
        display = data[columns].rename(
            columns={
                "scenario": "Scenario",
                "return_period_years": "Return period (yr)",
                "R_mm": "Rainfall (mm)",
                "H_m": f"Water level (m, {settings.datum})",
                "H_model_m": "Model water level (m)",
                "diagnostic_only": "Diagnostic only",
            }
        )
        table = display.to_html(
            index=False, border=0, float_format=lambda x: f"{x:.3f}", na_rep="n.a.", escape=True
        )
        used = [r for r in sources if r["unit"] == meta["unit"] and r["selected_observations"] > 0]
        source_text = ", ".join(escape(r["station_id"]) for r in used)
        sections.append(f"""<section id="unit-{escape(meta["unit"])}"><div class="eyebrow">Analysis unit</div><h2>{escape(meta["unit"])} · {escape(meta["gauge"])}</h2>
<p class="muted">{meta["event_count"]} events · {meta["effective_years"]:.2f} effective paired years · {escape(meta["copula_family"])} copula · Joint CvM p = {meta["CvM_p"]:.3f}</p>
<figure><img src="units/{escape(meta["unit"])}/figures/source_coverage.svg" alt="Monthly coverage of contributing rainfall observations"><figcaption>Only rainfall sources with a real contribution are shown. Empty coverage is not interpreted as zero rainfall.</figcaption></figure>
<figure><img src="units/{escape(meta["unit"])}/figures/model_diagnostics.svg" alt="Marginal distributions and Kendall design curves"><figcaption>Stars identify maximum physical joint-density design points. Curves describe events selected using the rainfall threshold.</figcaption></figure>
<figure><img src="units/{escape(meta["unit"])}/figures/annual_frequency.svg" alt="Full-record rainfall and tide-station annual frequency curves"><figcaption>S1 and S2 use full, separately screened records. S2 is fitted once per tide station. Gray points are accepted annual observations; stars identify design levels.</figcaption></figure>
<div class="scroll">{table}</div><p><small>Contributing rainfall stations: {source_text}</small></p></section>""")
    notice = (
        '<p class="notice"><strong>Synthetic demonstration.</strong> All observations are generated. These values are not study estimates.</p>'
        if settings.synthetic
        else ""
    )
    nav = "".join(f'<a href="#unit-{escape(m["unit"])}">{escape(m["unit"])}</a>' for m in models)
    body = f"""<header><div class="eyebrow">Joint return-period analysis</div><h1>{title}</h1><p class="muted">Compound return periods and full-record univariate designs</p></header>{notice}
<div class="grid"><div class="stat"><b>{len(models)}</b>compound models</div><div class="stat"><b>{sum(not m["diagnostic_only"] for m in models)}</b>models meeting formal gates</div><div class="stat"><b>{settings.datum}</b>analysis vertical reference</div></div>
<section><h2>Design definitions</h2><p>S1 uses full-record annual rainfall maxima. S2 uses full-record annual coastal water-level maxima, fitted once per tide station. Compound designs use ordinary Kendall return periods and the maximum physical joint-density point.</p>
<p><a href="design_values.csv">Design values</a> · <a href="rainfall_sources_used.csv">Contributing rainfall stations</a> · <a href="provenance.json">Input provenance</a> · <a href="settings.json">Numerical settings</a></p><nav>{nav}</nav></section>{"".join(sections)}
<footer>Point estimates. No parameter-uncertainty intervals are supplied. Diagnostic status is retained explicitly.</footer>"""
    html = f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title} | Return-period analysis</title><style>{CSS}</style></head><body><main>{body}</main></body></html>'
    (out / "results.html").write_text(html, encoding="utf8")


def annual_figure(out, meta, settings):
    style()
    folder = out / "units" / meta["unit"]
    rain = json.loads((folder / "rain_annual_model.json").read_text(encoding="utf8"))
    tides = json.loads((out / "tide_models.json").read_text(encoding="utf8"))
    tide = tides[meta["gauge"]]
    transfer = tide.get("transfer", {})
    donor = tides.get(transfer.get("donor"), tide)
    designs = pd.read_csv(out / "design_values.csv", dtype={"unit": str})
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.3), layout="constrained")
    periods = np.geomspace(2, 200, 240)
    for j, (variable, table, model) in enumerate(
        [
            ("R", folder / "rain_annual_maxima.csv", rain),
            ("H", out / "tides" / meta["gauge"] / "annual_maxima.csv", donor),
        ]
    ):
        ax = axes[j]
        annual = pd.read_csv(table)
        values = np.sort(annual.loc[annual.used, "maximum"].to_numpy(float))
        factor = settings.rainfall_factor_24h if variable == "R" else 1.0
        if len(values):
            empirical = 1 / (1 - (np.arange(len(values)) + 0.5) / len(values))
            ax.scatter(
                empirical,
                values * factor,
                s=15,
                color=COLORS["sample"],
                label="Accepted annual sample",
            )
        curve = None
        if "spec" in model:
            curve = Margin(model["spec"]).ppf(1 - 1 / periods) * factor
            if variable == "H":
                curve += transfer.get("delta_m", 0.0)
        elif variable == "R" and "L_CV" in model:
            from .regional import growth_curve

            growth, _ = growth_curve(model["L_CV"], model["L_skew"], periods)
            curve = np.array(list(growth.values())) * model["index_mm"] * factor
        if curve is not None:
            ax.plot(
                periods,
                curve,
                color=COLORS["rain" if variable == "R" else "water"],
                label="Annual-frequency curve",
            )
        selected = designs.loc[designs.scenario.eq("S1" if variable == "R" else "S2")]
        selected = (
            selected.loc[selected.unit.eq(meta["unit"])]
            if variable == "R"
            else selected.loc[selected.gauge.eq(meta["gauge"])]
        )
        ax.scatter(
            selected.return_period_years,
            selected.R_mm if variable == "R" else selected.H_m,
            marker="*",
            s=90,
            color=COLORS["rain" if variable == "R" else "water"],
            edgecolors="white",
            linewidths=0.4,
            zorder=4,
        )
        ax.set(
            xscale="log",
            xlabel="Return period (yr)",
            ylabel="24 h rainfall (mm)"
            if variable == "R"
            else f"Coastal water level (m, {settings.datum})",
        )
        ax.set_title(
            "a  S1 rainfall" if variable == "R" else "b  S2 coastal water level", loc="left"
        )
        ax.set_xticks([2, 10, 50, 100, 200], ["2", "10", "50", "100", "200"])
        ax.grid(alpha=0.8)
        ax.legend(fontsize=7)
    save(fig, folder / "figures/annual_frequency")
