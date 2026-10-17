#!/usr/bin/env python3
"""Reproduce Figure 4 from its adjacent frozen, revised source table.

Usage: python plot.py
Requires Python 3 and Matplotlib. No enrichment analyses are rerun.
The 23 frozen rows omit redundant KEGG map05204. The PDF is 180 x 210 mm.
"""

from pathlib import Path
import csv
import math
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


import argparse
HERE = Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--outdir", type=Path, default=HERE)
args=parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)
SOURCE = HERE / "Figure_4_source_data.tsv"
STEM = args.outdir / "Figure_4"
PANEL_CONFIG = {
    "A": ("GP-up GO terms", "#0072B2", 12.5, 2, 8),
    "B": ("GP-up KEGG pathways", "#0072B2", 32.5, 5, 4),
    "C": ("UP-up GO terms", "#D55E00", 31.0, 5, 8),
    "D": ("UP-up KEGG pathways", "#D55E00", 14.5, 2, 3),
}


def load_data():
    panels = {key: [] for key in PANEL_CONFIG}
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["gene_count"] = int(row["gene_count"])
            row["eligible_genes"] = int(row["eligible_genes"])
            row["gene_ratio_percent"] = float(row["gene_ratio_percent"])
            row["order"] = int(row["order"])
            if row["term_id"] == "map05204":
                raise ValueError("The redundant map05204 row is not in this figure.")
            if not math.isclose(row["gene_ratio_percent"],
                                100 * row["gene_count"] / row["eligible_genes"],
                                abs_tol=1e-8):
                raise ValueError("Gene ratio is inconsistent with its source counts.")
            panels[row["panel"]].append(row)
    for key, rows in panels.items():
        if len(rows) != PANEL_CONFIG[key][4]:
            raise ValueError(f"Unexpected frozen row count in panel {key}.")
        rows.sort(key=lambda row: row["order"])
    return panels


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.unicode_minus": False,
    })
    panels = load_data()
    width_mm, height_mm = 180, 210
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    top_mm = 12
    for panel, (title, color, xmax, major, _) in PANEL_CONFIG.items():
        rows = panels[panel]
        plot_height_mm = len(rows) * 6
        ax = fig.add_axes([
            80 / width_mm,
            1 - (top_mm + plot_height_mm) / height_mm,
            93 / width_mm,
            plot_height_mm / height_mm,
        ])
        positions = list(range(len(rows)))
        values = [row["gene_ratio_percent"] for row in rows]
        labels = [textwrap.fill(row["label"], width=43,
                               break_long_words=False, break_on_hyphens=False)
                  for row in rows]
        ax.barh(positions, values, height=0.62, color=color, zorder=3)
        ax.set_yticks(positions, labels=labels)
        for label in ax.get_yticklabels():
            label.set_linespacing(1.0)
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_xlim(0, xmax)
        ax.xaxis.set_major_locator(MultipleLocator(major))
        ax.tick_params(axis="y", labelsize=8.5, length=0, pad=5)
        ax.tick_params(axis="x", labelsize=8.5, length=3, color="#666666")
        ax.grid(axis="x", color="#DDDDDD", linewidth=0.5, zorder=0)
        for name in ("top", "right", "left"):
            ax.spines[name].set_visible(False)
        ax.spines["bottom"].set_color("#AAAAAA")
        ax.spines["bottom"].set_linewidth(0.6)
        for y, row in enumerate(rows):
            ax.text(row["gene_ratio_percent"] + xmax * 0.012, y,
                    f"n={row['gene_count']}", va="center", ha="left", fontsize=8)
        fig.text(5 / width_mm, 1 - (top_mm - 4) / height_mm,
                 panel.lower(), fontsize=11, fontweight="bold", va="baseline")
        fig.text(13 / width_mm, 1 - (top_mm - 4) / height_mm,
                 title, fontsize=10.5, fontweight="bold", va="baseline")
        top_mm += plot_height_mm + 14
    fig.text(126.5 / width_mm, 5 / height_mm, "Gene ratio (%)",
             fontsize=9, ha="center", va="bottom")
    fig.savefig(STEM.with_suffix(".pdf"), facecolor="white",
                metadata={"Title": "Figure 4. Plant functional enrichment"})
    fig.savefig(STEM.with_suffix(".png"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Reproduced 23 frozen rows from {SOURCE.name}; no analyses rerun.")


if __name__ == "__main__":
    main()
