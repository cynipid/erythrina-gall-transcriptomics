#!/usr/bin/env python3
"""Rebuild this figure using only the files beside this script.

Run: python plot.py
Dependencies: Python 3, matplotlib, numpy.
This is a rendering script: no DESeq2 fitting, PCA recomputation, shrinkage,
enrichment, or significance-threshold changes are performed.
"""
from pathlib import Path
import csv
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import argparse

HERE = Path(__file__).resolve().parent
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--outdir', type=Path, default=None, help='Output directory; source inputs remain beside this script.')
_args = _parser.parse_args()
OUTDIR = _args.outdir.expanduser().resolve() if _args.outdir is not None else HERE
OUTDIR.mkdir(parents=True, exist_ok=True)
META = json.loads((HERE / "metadata.json").read_text())
with (HERE / META["source_data"]).open(newline="") as handle:
    ROWS = list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    """Explicit float parsing preserves very small fixed-point P values."""
    return float(value) if value not in ("", "NA", "NaN", "nan") else math.nan


def setup():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.fontsize": 10, "axes.linewidth": 0.9,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(180 / 25.4, META["height_mm"] / 25.4))
    fig.subplots_adjust(left=0.13, right=0.975, bottom=0.13, top=0.86)
    for edge in ["top", "right"]:
        ax.spines[edge].set_visible(False)
    ax.tick_params(length=4, width=0.9)
    return fig, ax


def draw_pca(fig, ax):
    assert len(ROWS) == 6
    handles = []
    for group in META["groups"]:
        selected = [r for r in ROWS if r["group"] == group["key"]]
        assert len(selected) == 3
        x = [number(r["PC1"]) for r in selected]
        y = [number(r["PC2"]) for r in selected]
        ax.scatter(x, y, s=58, c=group["color"], marker=group["marker"],
                   linewidths=0.7, edgecolors="white", zorder=4)
        handles.append(Line2D([], [], color="none", marker=group["marker"],
                              markerfacecolor=group["color"], markeredgecolor="none",
                              markersize=7, label=group["label"]))
    for row in ROWS:
        dx, dy = number(row["label_dx_points"]), number(row["label_dy_points"])
        ax.annotate(row["sample_label"], (number(row["PC1"]), number(row["PC2"])),
                    xytext=(dx, dy), textcoords="offset points", ha="center",
                    va="center", fontsize=10, color="#222222",
                    arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "#777777",
                                "shrinkA": 2, "shrinkB": 6} if row["leader_line"] == "TRUE" else None,
                    zorder=5)
    ax.set_xlabel(f'PC1: {META["PC1_variance_percent"]:.2f}% variance')
    ax.set_ylabel(f'PC2: {META["PC2_variance_percent"]:.2f}% variance')
    ax.set_xlim(META["x_limits"])
    ax.set_ylim(META["y_limits"])
    ax.set_xticks(META["x_ticks"])
    ax.set_yticks(META["y_ticks"])
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.07),
              ncol=2, frameon=False, handletextpad=0.35, columnspacing=1.8)


def draw_volcano(fig, ax):
    counts = {g["key"]: sum(r["status"] == g["key"] for r in ROWS) for g in META["groups"]}
    assert counts == META["expected_status_counts"], (counts, META["expected_status_counts"])
    for row in ROWS:
        lfc, padj = number(row["log2FoldChange"]), number(row["padj"])
        expected = "Not tested" if not math.isfinite(padj) else "Non-DE"
        if math.isfinite(padj) and padj < 0.01 and abs(lfc) >= 1:
            expected = META["positive_status"] if lfc > 0 else META["negative_status"]
        assert row["status"] == expected, row["gene_id"]
    handles = []
    for group in META["groups"]:
        selected = [r for r in ROWS if r["status"] == group["key"] and r["plotted"] == "TRUE"]
        if not selected:
            continue
        x = [number(r["log2FoldChange"]) for r in selected]
        y = [number(r["plot_y"]) for r in selected]
        marker = group["marker"]
        ax.scatter(x, y, s=8 if marker == "x" else 5, color=group["color"], marker=marker,
                   alpha=0.8 if marker == "x" else 0.62,
                   linewidths=0.6 if marker == "x" else 0,
                   zorder=group["zorder"])
        handles.append(Line2D([], [], linestyle="none", marker=marker,
                              color=group["color"], markersize=5,
                              label=f'{group["label"]} (n = {counts[group["key"]]:,})'))
    for x in [-1, 1]:
        ax.axvline(x, color="#777777", lw=0.8, ls=(0, (4, 4)), zorder=0)
    ax.axhline(2, color="#777777", lw=0.8, ls=(0, (4, 4)), zorder=0)
    ax.set_xlabel(META["x_label"])
    ax.set_ylabel(r'$-\log_{10}$ adjusted P')
    ax.set_xlim(META["x_limits"])
    ax.set_ylim(META["y_limits"])
    ax.set_xticks(META["x_ticks"])
    ax.set_yticks(META["y_ticks"], labels=META["y_tick_labels"])
    order = META["legend_order"]
    handles.sort(key=lambda h: order.index(h.get_label().split(" (n =")[0]))
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.025),
              ncol=2, frameon=False, handletextpad=0.4, columnspacing=1.7)


def main():
    fig, ax = setup()
    if META["kind"] == "PCA":
        draw_pca(fig, ax)
    else:
        draw_volcano(fig, ax)
    stem = META["output_stem"]
    fig.savefig(OUTDIR / f"{stem}.pdf", metadata={"Creator": "Matplotlib",
                "CreationDate": None, "ModDate": None})
    fig.savefig(OUTDIR / f"{stem}.png", dpi=300)
    plt.close(fig)
    print(f'{stem}: reproduced from {len(ROWS):,} source rows.')


if __name__ == "__main__":
    main()
