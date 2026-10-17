#!/usr/bin/env python3
"""Reproduce Figure 2 from its adjacent frozen source table.

Usage: python plot.py
Requires Python 3 and Matplotlib. No enrichment analyses are rerun.
The PDF has a fixed 180 x 240 mm page; fonts are sized for that print width.
"""

from pathlib import Path
import csv

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
SOURCE = HERE / "Figure_2_source_data.tsv"
STEM = args.outdir / "Figure_2"
PANEL_CONFIG = {
    "A": ("Larva-up GO terms", "#009E73"),
    "B": ("Larva-up KEGG pathways", "#009E73"),
    "C": ("Pupa-up GO terms", "#CC79A7"),
    "D": ("Pupa-up KEGG pathways", "#CC79A7"),
}


def load_data():
    panels = {key: [] for key in PANEL_CONFIG}
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["overlap"] = int(row["overlap"])
            row["gene_ratio_percent"] = float(row["gene_ratio_percent"])
            row["padj"] = float(row["padj"])
            panels[row["panel"]].append(row)
    if any(len(rows) != 8 for rows in panels.values()):
        raise ValueError("Expected the frozen eight rows in each panel.")
    # The earlier plot sorted ascending on a non-inverted y axis. Reversing
    # that stable order preserves the exact top-to-bottom order, including ties.
    return {
        key: sorted(rows, key=lambda row: row["gene_ratio_percent"])[::-1]
        for key, rows in panels.items()
    }


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.unicode_minus": False,
    })
    panels = load_data()
    width_mm, height_mm = 180, 240
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    for index, (panel, (title, color)) in enumerate(PANEL_CONFIG.items()):
        # Full-width panels make complete term names readable at 180 mm.
        top_mm = 12 + index * 57
        plot_height_mm = 40
        ax = fig.add_axes([
            80 / width_mm,
            1 - (top_mm + plot_height_mm) / height_mm,
            93 / width_mm,
            plot_height_mm / height_mm,
        ])
        rows = panels[panel]
        positions = list(range(len(rows)))
        values = [row["gene_ratio_percent"] for row in rows]
        ax.barh(positions, values, height=0.63, color=color, zorder=3)
        ax.set_yticks(positions, labels=[row["term_name"] for row in rows])
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_xlim(0, 16)
        ax.xaxis.set_major_locator(MultipleLocator(4))
        ax.tick_params(axis="y", labelsize=8.5, length=0, pad=5)
        ax.tick_params(axis="x", labelsize=8.5, length=3, color="#666666")
        ax.grid(axis="x", color="#DDDDDD", linewidth=0.5, zorder=0)
        for name in ("top", "right", "left"):
            ax.spines[name].set_visible(False)
        ax.spines["bottom"].set_color("#AAAAAA")
        ax.spines["bottom"].set_linewidth(0.6)
        for y, row in enumerate(rows):
            ax.text(row["gene_ratio_percent"] + 0.17, y,
                    f"n={row['overlap']}", va="center", ha="left", fontsize=8)
        fig.text(5 / width_mm, 1 - (top_mm - 4) / height_mm,
                 panel.lower(), fontsize=11, fontweight="bold", va="baseline")
        fig.text(13 / width_mm, 1 - (top_mm - 4) / height_mm,
                 title, fontsize=10.5, fontweight="bold", va="baseline")
    fig.text(126.5 / width_mm, 5 / height_mm, "Gene ratio (%)",
             fontsize=9, ha="center", va="bottom")
    # No tight cropping: preserve the specified physical size and font sizes.
    fig.savefig(STEM.with_suffix(".pdf"), facecolor="white",
                metadata={"Title": "Figure 2. Insect functional enrichment"})
    fig.savefig(STEM.with_suffix(".png"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Reproduced 32 frozen rows from {SOURCE.name}; no analyses rerun.")


if __name__ == "__main__":
    main()
