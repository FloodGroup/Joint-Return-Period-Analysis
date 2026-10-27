"""Shared publication figure style for every study region."""

from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = {
    "ink": "#203844",
    "muted": "#71818a",
    "grid": "#e5eaed",
    "rain": "#168c98",
    "water": "#4c75a3",
    "diagnostic": "#ba7143",
    "sample": "#b9c4ca",
}
PERIOD_COLORS = ("#168c98", "#4c75a3", "#8b68a1", "#ba7143")


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelcolor": COLORS["ink"],
            "text.color": COLORS["ink"],
            "xtick.color": COLORS["muted"],
            "ytick.color": COLORS["muted"],
            "axes.edgecolor": "#aebdc4",
            "axes.spines.top": True,
            "axes.spines.right": True,
            "axes.spines.bottom": True,
            "axes.spines.left": True,
            "axes.linewidth": 0.8,
            "grid.color": COLORS["grid"],
            "grid.linewidth": 0.65,
            "axes.axisbelow": True,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "lines.linewidth": 1.8,
            "figure.dpi": 140,
            "savefig.dpi": 300,
        }
    )


def save(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".svg", ".png"):
        fig.savefig(path.with_suffix(suffix), bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
