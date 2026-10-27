"""04. Fit event margins and copula dependence."""

import re

from jrpa.common.io import write_csv, write_json

from .dependence import fit_dependence


def step04_copula(analysis, skip_compound=False):
    """04. Fit event margins and copulas; retain goodness-of-fit diagnostics."""
    analysis.require(4)
    settings = analysis.settings
    for item in analysis.units:
        uid, folder = item["unit"]["unit"], item["folder"]
        pot, event_meta = item["pot"], item["event_meta"]
        if not skip_compound:
            number = int(re.sub(r"\D", "", uid) or 0)
            meta, logs = fit_dependence(
                pot,
                event_meta,
                settings,
                number,
                water_override=settings.region == "gba" and uid == "10a" and not settings.synthetic,
            )
            write_json(folder / "compound_model.json", meta)
            for label, rows in logs.items():
                write_csv(folder / (label + ".csv"), rows)
            item["compound_model"] = meta
            analysis.fitted.append(meta)
    analysis.completed = 4
