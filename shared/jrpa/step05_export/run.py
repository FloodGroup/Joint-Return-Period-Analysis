"""05. Calculate Kendall design points and export tables and figures."""

import re

from jrpa.common.datums import model_offset
from jrpa.common.io import sha256, write_csv, write_json

from .kendall import design_points
from .reporting import render_report


def step05_export(analysis):
    """05. Calculate Kendall design points and export tables and figures."""
    analysis.require(5)
    settings, out = analysis.settings, analysis.out
    design_rows, fitted, sources = analysis.design_rows, analysis.fitted, analysis.sources
    annual_tide, inputs, tide_quality = analysis.annual_tide, analysis.inputs, analysis.tide_quality
    manifest_path = analysis.manifest_path
    for item in analysis.units:
        uid, gauge = item["unit"]["unit"], item["unit"]["gauge"]
        folder = item["folder"]
        if "compound_model" in item:
            meta = item["compound_model"]
            number = int(re.sub(r"\D", "", uid) or 0)
            points, curves, calibration = design_points(meta, settings, number)
            for point in points:
                point["gauge"] = gauge
                point["H_model_m"] = point["H_m"] + model_offset(settings.region, gauge)
            write_csv(folder / "kendall_design.csv", points)
            write_csv(folder / "kendall_curves.csv", curves)
            write_json(folder / "kendall_calibration.json", calibration)
            design_rows.extend(points)
    analysis.completed = 5
    write_csv(out / "design_values.csv", design_rows)
    write_csv(out / "rainfall_sources_used.csv", sources)
    write_json(out / "tide_models.json", annual_tide)
    write_json(out / "settings.json", settings.as_dict())
    write_json(
        out / "provenance.json",
        dict(
            manifest_sha256=sha256(manifest_path),
            inputs=inputs,
            tide_quality=tide_quality,
            synthetic=settings.synthetic,
        ),
    )
    render_report(out, settings, design_rows, fitted, sources)
    return dict(
        region=settings.region,
        units=len(analysis.units),
        design_records=len(design_rows),
        output=str(out),
        synthetic=settings.synthetic,
    )
