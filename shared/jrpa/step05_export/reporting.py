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
:root{--ink:#24343d;--muted:#566872;--accent:#087f86;--line:#dce4e8;--paper:#fff}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:88px}
body{margin:0;color:var(--ink);font:16px/1.65 "Segoe UI",Arial,sans-serif;background:#f4f6f7}
main{max-width:1240px;margin:auto;padding:40px 32px 80px}
header{padding:24px 0 28px;border-top:3px solid var(--accent)}
h1{font-size:clamp(30px,4vw,44px);line-height:1.15;letter-spacing:-.025em;margin:10px 0 14px}
h2{font-size:28px;line-height:1.25;margin:4px 0 12px}h3{font-size:20px;line-height:1.4;margin:32px 0 16px}
p{margin:10px 0 16px}.eyebrow{font-size:12px;font-weight:700;letter-spacing:.12em;color:var(--accent);text-transform:uppercase}
.muted,small{color:var(--muted)}small{font-size:14px}a{color:#076c75;text-underline-offset:3px}
a:hover{color:#034d54}a:focus-visible,summary:focus-visible,.scroll:focus-visible{outline:3px solid #cf8434;outline-offset:4px}
.unit,.overview{background:var(--paper);border:1px solid var(--line);padding:32px;margin:28px 0;border-radius:8px}
.unit-head{display:flex;justify-content:space-between;align-items:baseline;gap:16px;flex-wrap:wrap;border-bottom:1px solid var(--line);padding-bottom:16px}
.unit-head p{margin:0}.notice{background:#fff4e6;border-left:3px solid #ad6c25;padding:12px 18px}
.stats{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:20px;margin:8px 0 28px}
.stat{border-top:1px solid var(--line);padding:14px 0}.stat b{display:block;font-size:24px;line-height:1.4}.stat span{font-size:14px;color:var(--muted)}
nav{display:flex;gap:8px;flex-wrap:wrap;position:sticky;top:0;z-index:10;background:#f4f6f7f5;padding:12px 0;border-bottom:1px solid var(--line)}
nav a{padding:6px 12px;background:white;border:1px solid #cbd8dd;border-radius:4px;text-decoration:none;font-size:14px}
.downloads{display:flex;gap:10px 24px;flex-wrap:wrap;font-size:14px;margin-top:20px}
.figure-grid,.annual-tables{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:24px;align-items:start}
.annual-tables>*{min-width:0}.annual-tables th,.annual-tables td{padding:12px}
figure{margin:0;min-width:0}figure img{display:block;width:100%;height:auto}figure a{display:block}
figcaption{font-size:14px;color:var(--muted);line-height:1.6;margin:10px 4px 0}
.kendall{max-width:760px;margin:24px auto 0}.coverage{margin:20px 0}.coverage a{overflow-x:auto}.coverage img{min-width:640px}
.scroll{overflow-x:auto;border:1px solid var(--line);border-radius:5px;margin:12px 0 20px}
table{border-collapse:collapse;width:100%;font-size:15px;line-height:1.5;font-variant-numeric:tabular-nums}
caption{text-align:left;font-weight:600;font-size:16px;padding:14px 18px;background:#f0f5f6;caption-side:top}
th,td{padding:12px 18px;border-bottom:1px solid var(--line);text-align:right;vertical-align:top}
th{font-size:14px;font-weight:600;background:#f8fafb;color:#425762}th:first-child,td:first-child{text-align:left}
thead th{vertical-align:bottom}tbody tr:last-child td{border-bottom:0}tbody tr:nth-child(even){background:#fafcfc}
tbody tr:hover{background:#edf5f6}td{white-space:nowrap}.unit-symbol{display:block;font-weight:400;font-size:13px}
.status{color:#8a541b}.details-table td:last-child{white-space:normal;min-width:240px;text-align:left}
details{border-top:1px solid var(--line);margin-top:24px;padding-top:16px}summary{cursor:pointer;font-size:15px;font-weight:600}
.fit-meta{display:flex;gap:10px 28px;flex-wrap:wrap;margin:12px 0 20px;color:var(--muted);font-size:14px}
.fit-meta b{color:var(--ink)}.empty{padding:18px;background:#f8fafb;color:var(--muted)}
@media(max-width:980px){.figure-grid,.annual-tables{grid-template-columns:1fr}}
@media(max-width:760px){main{padding:20px 16px 40px}.unit,.overview{padding:20px 16px;margin:20px 0}.stats{grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}h2{font-size:24px}h3{font-size:19px}nav{position:static}th,td{padding:10px 12px}.unit-head{gap:4px}table{min-width:400px}}
@media print{body{background:white}main{max-width:none;padding:0}nav,.downloads{display:none}.unit,.overview{padding:20px 0;border:0;border-radius:0}.unit{break-before:page}.scroll{overflow:visible}figure,table{break-inside:avoid}.figure-grid{grid-template-columns:repeat(2,minmax(0,1fr))}a{color:inherit;text-decoration:none}}
"""


def unit_figures(root, meta, synthetic=False):
    style()
    plt.rcParams.update(
        {"font.size": 12, "axes.labelsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 11}
    )
    root = Path(root)
    data = pd.read_csv(root / "pot_events.csv")
    curves = pd.read_csv(root / "kendall_curves.csv")
    points = pd.read_csv(root / "kendall_design.csv")
    for j, (variable, col, label) in enumerate(
        [
            ("R", "rain_mm", "Event rainfall (mm)"),
            ("H", "water_max_m", f"Coastal water level (m, {meta['datum']})"),
        ]
    ):
        x = np.sort(data[col].to_numpy(float))
        f = (np.arange(len(x)) + 0.5) / len(x)
        fig, ax = plt.subplots(figsize=(5.8, 4.4), layout="constrained")
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
        ax.set_title("a  Rainfall margin" if j == 0 else "b  Water-level margin", loc="left")
        ax.grid(alpha=0.8)
        ax.legend(fontsize=10, loc="lower right")
        save(fig, root / f"figures/margin_{variable}")
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
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
    ax.legend(fontsize=10)
    save(fig, root / "figures/kendall_design")


def report_figure(uid, name, label, caption="", css=""):
    """Link each vector figure to its full-size version."""
    src = f"units/{escape(uid)}/figures/{name}.svg"
    return (
        f'<figure class="{css}"><a href="{src}" target="_blank" rel="noopener" '
        f'aria-label="Open {escape(label)} at full size">'
        f'<img src="{src}" alt="{escape(label)}" loading="lazy"></a>'
        f"<figcaption>{escape(caption)}</figcaption></figure>"
    )


def design_table(data, scenario, datum):
    """Show each design class with its relevant variables and units."""
    part = data.loc[data.scenario.eq(scenario)].sort_values("return_period_years")
    if part.empty:
        return ""
    labels = {
        "S1": "S1 · Univariate design rainfall",
        "S2": "S2 · Univariate design coastal water level",
        "Kendall": "Compound · Kendall",
    }
    columns = [("return_period_years", "Return period", "yr", 0)]
    if scenario != "S2":
        columns.append(("R_mm", "24 h rainfall", "mm", 2))
    if scenario != "S1":
        columns.extend(
            [
                ("H_m", "Coastal water level", f"m, {datum}", 3),
                ("H_model_m", "Model water level", "m", 3),
            ]
        )
    diagnostic = part.get("diagnostic_only", pd.Series(False, index=part.index)).fillna(False)
    show_status = diagnostic.any()
    heads = "".join(
        f'<th scope="col">{escape(label)}<span class="unit-symbol">{escape(unit)}</span></th>'
        for _, label, unit, _ in columns
    )
    if show_status:
        heads += '<th scope="col">Status</th>'
    rows = []
    for index, row in part.iterrows():
        cells = "".join(
            f"<td>{float(row[key]):,.{digits}f}</td>" if pd.notna(row.get(key)) else "<td>—</td>"
            for key, _, _, digits in columns
        )
        if show_status:
            cells += (
                '<td class="status">Diagnostic</td>'
                if diagnostic.loc[index]
                else "<td>Accepted</td>"
            )
        rows.append(f"<tr>{cells}</tr>")
    label = labels.get(scenario, scenario)
    return (
        f'<div class="scroll" tabindex="0" role="region" aria-label="{escape(label)} design values">'
        f"<table><caption>{escape(label)}</caption><thead><tr>{heads}</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def render_report(output, settings, design_rows, models, sources, units=None):
    out = Path(output)
    title = {"gba": "Greater Bay Area", "san_francisco": "San Francisco Bay", "tokyo": "Tokyo Bay"}[
        settings.region
    ]
    frame = pd.DataFrame(design_rows).reindex(
        columns=[
            "unit",
            "gauge",
            "scenario",
            "return_period_years",
            "R_mm",
            "H_m",
            "H_model_m",
            "diagnostic_only",
            "diagnostic_reason",
        ]
    )
    model_by_unit = {m["unit"]: m for m in models}
    if units is None:
        # Support regeneration from previously exported tables and model files.
        units = (
            frame.loc[frame.unit.notna(), ["unit", "gauge"]].drop_duplicates().to_dict("records")
        )
        known = {u["unit"] for u in units}
        units.extend(
            {"unit": m["unit"], "gauge": m["gauge"]} for m in models if m["unit"] not in known
        )
    sections = []
    for unit in units:
        uid, gauge = unit["unit"], unit["gauge"]
        folder = out / "units" / uid
        meta = model_by_unit.get(uid)
        annual_figure(out, unit, settings)
        data = frame.loc[
            frame.unit.eq(uid) | (frame.scenario.eq("S2") & frame.gauge.eq(gauge))
        ].copy()
        annual_tables = "".join(
            design_table(data, scenario, settings.datum) for scenario in ("S1", "S2")
        )
        tables = '<div class="annual-tables">' + annual_tables + "</div>" if annual_tables else ""
        tables += design_table(data, "Kendall", settings.datum)
        if not tables:
            tables = '<p class="empty">No design values were produced for this unit.</p>'
        annual = (
            '<div class="figure-grid">'
            + report_figure(
                uid,
                "annual_R",
                "S1 rainfall frequency curve",
                "S1 · Annual maxima of rainfall and estimated 24-hour design rainfall.",
            )
            + report_figure(
                uid,
                "annual_H",
                "S2 coastal water-level frequency curve",
                "S2 · Annual maxima of coastal water level and estimated design water levels.",
            )
            + "</div>"
        )
        compound = ""
        if meta is not None:
            unit_figures(folder, meta, settings.synthetic)
            compound = (
                '<h3>Compound analysis</h3><div class="fit-meta">'
                f"<span><b>{meta['event_count']}</b> events</span>"
                f"<span><b>{meta['effective_years']:.2f}</b> paired years</span>"
                f"<span><b>{escape(meta['copula_family'])}</b> copula</span>"
                f"<span>Joint CvM <i>p</i> = <b>{meta['CvM_p']:.3f}</b></span></div>"
                '<div class="figure-grid">'
                + report_figure(uid, "margin_R", "Event rainfall distribution")
                + report_figure(uid, "margin_H", "Event water-level distribution")
                + "</div>"
                + report_figure(
                    uid,
                    "kendall_design",
                    "Kendall return-period contours",
                    "Contours show Kendall return periods; stars mark maximum joint-density design combinations.",
                    "kendall",
                )
            )
        used = [r for r in sources if r["unit"] == uid and r["selected_observations"] > 0]
        source_text = ", ".join(escape(r["station_id"]) for r in used)
        coverage = ""
        if (folder / "figures/source_coverage.svg").exists():
            coverage = "<h3>Observation coverage</h3>" + report_figure(
                uid,
                "source_coverage",
                "Coverage of contributing rainfall observations",
                "Monthly fraction of observations supplied by each rainfall station.",
                "coverage",
            )
        diagnostic = data.loc[
            data.diagnostic_only.fillna(False).astype(bool), ["scenario", "diagnostic_reason"]
        ].drop_duplicates()
        details = ""
        if not diagnostic.empty:
            details = (
                '<details><summary>Fit diagnostics</summary><div class="scroll details-table">'
                + diagnostic.rename(
                    columns={"scenario": "Design", "diagnostic_reason": "Reason"}
                ).to_html(index=False, border=0, escape=True, na_rep="—")
                + "</div></details>"
            )
        sections.append(
            f'<section class="unit" id="unit-{escape(uid)}"><div class="unit-head"><div>'
            f'<div class="eyebrow">Analysis unit</div><h2>Unit {escape(uid)}</h2></div>'
            f'<p class="muted">Tide station <b>{escape(gauge)}</b></p></div>'
            f"<h3>Design values</h3>{tables}<h3>Annual frequency analysis</h3>{annual}"
            '<p class="muted"><small>Rainfall and water level use separately screened annual maxima. '
            "Points show annual observations; stars show design values.</small></p>"
            f"{compound}{coverage}<p><small>Contributing rainfall stations: {source_text or 'None'}</small></p>"
            f"{details}</section>"
        )
    notice = (
        '<p class="notice"><strong>Synthetic demonstration.</strong> Generated observations.</p>'
        if settings.synthetic
        else ""
    )
    nav = "".join(
        f'<a href="#unit-{escape(u["unit"])}">Unit {escape(u["unit"])}</a>' for u in units
    )
    periods = ", ".join(str(t) for t in settings.periods)
    mode = (
        "Univariate and compound return-period analysis"
        if models
        else "Univariate return-period analysis"
    )
    summary = [
        (str(len(units)), "Analysis units"),
        (str(len({u["gauge"] for u in units})), "Tide stations"),
        (periods, "Return periods (yr)"),
        (settings.datum, "Vertical reference"),
    ]
    stats = "".join(
        f'<div class="stat"><b>{escape(value)}</b><span>{label}</span></div>'
        for value, label in summary
    )
    definition = "S1 denotes univariate 24-hour design rainfall; S2 denotes univariate design coastal water level. Both are estimated from separately screened annual-maximum samples."
    if models:
        definition += " Compound combinations are selected on Kendall return-period contours using the maximum joint-density criterion."
    body = (
        f'<header><div class="eyebrow">Return-period analysis</div><h1>{title}</h1><p class="muted">{mode}</p></header>'
        f'{notice}<div class="stats">{stats}</div><section class="overview"><h2>Design events</h2><p>{definition}</p>'
        '<div class="downloads"><a href="design_values.csv">Design values · CSV</a>'
        '<a href="rainfall_sources_used.csv">Rainfall stations · CSV</a>'
        '<a href="settings.json">Analysis settings</a><a href="provenance.json">Data sources</a></div></section>'
        f'<nav aria-label="Analysis units">{nav}</nav>{"".join(sections)}'
    )
    html = f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title} | Return-period analysis</title><style>{CSS}</style></head><body><main>{body}</main></body></html>'
    (out / "results.html").write_text(html, encoding="utf8")


def annual_figure(out, meta, settings):
    style()
    plt.rcParams.update(
        {"font.size": 12, "axes.labelsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 11}
    )
    folder = out / "units" / meta["unit"]
    rain = json.loads((folder / "rain_annual_model.json").read_text(encoding="utf8"))
    tides = json.loads((out / "tide_models.json").read_text(encoding="utf8"))
    tide = tides[meta["gauge"]]
    transfer = tide.get("transfer", {})
    donor = tides.get(transfer.get("donor"), tide)
    try:
        designs = pd.read_csv(out / "design_values.csv", dtype={"unit": str})
    except pd.errors.EmptyDataError:
        designs = pd.DataFrame()
    designs = designs.reindex(
        columns=["unit", "gauge", "scenario", "return_period_years", "R_mm", "H_m"]
    )
    periods = np.geomspace(2, max(200, max(settings.periods)), 240)
    for j, (variable, table, model) in enumerate(
        [
            ("R", folder / "rain_annual_maxima.csv", rain),
            ("H", out / "tides" / meta["gauge"] / "annual_maxima.csv", donor),
        ]
    ):
        fig, ax = plt.subplots(figsize=(5.8, 4.4), layout="constrained")
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
            from jrpa.step03_annual.regional import growth_curve

            growth, _ = growth_curve(model["L_CV"], model["L_skew"], periods)
            curve = np.array(list(growth.values())) * model["index_mm"] * factor
        if curve is not None:
            ax.plot(
                periods,
                curve,
                color=COLORS["rain" if variable == "R" else "water"],
                label="Univariate return-level curve",
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
            "a  Rainfall return levels (S1)" if variable == "R" else "b  Water-level return levels (S2)", loc="left"
        )
        ax.set_xticks([2, 10, 50, 100, 200], ["2", "10", "50", "100", "200"])
        ax.grid(alpha=0.8)
        ax.legend(fontsize=10, loc="upper left")
        save(fig, folder / f"figures/annual_{variable}")
