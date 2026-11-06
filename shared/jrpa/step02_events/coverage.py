"""Monthly fractions of observations that actually contribute rainfall."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from jrpa.common.plotting import COLORS, save, style


def coverage_figure(folder, clock, blocks=None, daily_sources=None):
    style()
    local = clock.tz_localize(None) if clock.tz is not None else clock
    months = local.to_period("M")
    all_months = pd.period_range(months.min(), months.max(), freq="M")
    rows = []
    labels = []
    if daily_sources is not None:
        stations = daily_sources.dropna().unique()
        for sid in stations:
            counts = (
                pd.Series(daily_sources.eq(sid).to_numpy(int), index=months).groupby(level=0).sum()
            )
            denominator = np.array([p.days_in_month for p in all_months])
            rows.append(counts.reindex(all_months, fill_value=0).to_numpy() / denominator)
            labels.append(str(sid))
    else:
        for sid, part in blocks.groupby("station_id", sort=False):
            occupied = np.zeros(len(clock), bool)
            for r in part.itertuples():
                occupied[r.start_pos : r.end_pos + 1] = True
            if not occupied.any():
                continue
            counts = pd.Series(occupied.astype(int), index=months).groupby(level=0).sum()
            denominator = np.array([p.days_in_month * 24 for p in all_months])
            rows.append(counts.reindex(all_months, fill_value=0).to_numpy() / denominator)
            labels.append(str(sid))
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(11.5, max(2.5, 0.45 * len(rows) + 1.2)), layout="constrained")
    cmap = LinearSegmentedColormap.from_list("observed_rain", ["#f2f6f7", COLORS["rain"]])
    im = ax.imshow(np.array(rows), vmin=0, vmax=1, aspect="auto", interpolation="none", cmap=cmap)
    ax.set_yticks(range(len(labels)), labels, fontsize=10)
    step = 2 if len(all_months) <= 144 else 5 if len(all_months) <= 480 else 10
    ticks = (
        [0]
        + [
            i
            for i, p in enumerate(all_months)
            if p.month == 1 and p.year % step == 0 and 12 <= i < len(all_months) - 12
        ]
        + [len(all_months) - 1]
    )
    ticks = sorted(set(ticks))
    date_format = "%Y-%m" if all_months[0].year == all_months[-1].year else "%Y"
    ax.set_xticks(ticks, [all_months[i].strftime(date_format) for i in ticks])
    ax.set_title("Contributing rainfall observations", loc="left")
    ax.set_xlabel("Calendar year (UTC for hourly observations)")
    fig.colorbar(im, ax=ax, shrink=0.8, pad=0.03, ticks=[0, 0.5, 1], label="Monthly fraction")
    save(fig, Path(folder) / "figures/source_coverage")
