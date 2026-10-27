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
    fig, ax = plt.subplots(figsize=(11.5, max(1.8, 0.35 * len(rows) + 0.8)), layout="constrained")
    cmap = LinearSegmentedColormap.from_list("observed_rain", ["#f2f6f7", COLORS["rain"]])
    im = ax.imshow(np.array(rows), vmin=0, vmax=1, aspect="auto", interpolation="none", cmap=cmap)
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ticks = [i for i, p in enumerate(all_months) if p.month == 1 and p.year % 5 == 0]
    if not ticks:
        ticks = [0, len(all_months) - 1]
    ax.set_xticks(ticks, [str(all_months[i].year) for i in ticks])
    ax.set_title("Contributing rainfall observations", loc="left")
    ax.set_xlabel("Calendar year (UTC for hourly observations)")
    fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02, label="Selected fraction of calendar month")
    save(fig, Path(folder) / "figures/source_coverage")
