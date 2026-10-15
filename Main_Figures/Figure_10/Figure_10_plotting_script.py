#!/usr/bin/env python3
"""Reproduce Figure 10 using only the adjacent frozen source TSV.

Requires Python 3, pandas, numpy, and matplotlib. Run from any working directory:
    python /path/to/Figure_10_plotting_script.py
No expression analysis, MEROPS search, or enrichment test is rerun.
"""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=HERE / "Figure_10_source_data.tsv")
parser.add_argument("--outdir", type=Path, default=HERE)
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)
source = pd.read_csv(args.source, sep="\t")
assert source.groupby("panel").size().to_dict() == {"A": 90, "B": 9, "C": 19, "D": 2}
A, B, C, D = [source[source.panel.eq(p)].copy() for p in "ABCD"]
assert A.category.nunique() == 18 and int(A["count"].sum()) == 149
assert B.gene_id.nunique() == 9 and C.gene_id.nunique() == 19
assert B.plotted_log2FoldChange.le(-1).all() and C.plotted_log2FoldChange.ge(1).all()
assert np.allclose(B.original_log2FoldChange, B.plotted_log2FoldChange)
assert np.allclose(B.legacy_display_log2FoldChange, -B.plotted_log2FoldChange)
assert B.original_contrast.eq("pupa/larva").all() and B.plotted_contrast.eq("pupa/larva").all()
assert C.original_contrast.eq("GP/UP").all() and C.plotted_contrast.eq("GP/UP").all()
assert B.padj.lt(0.01).all() and C.padj.lt(0.01).all()
assert np.isfinite(B.lfcSE).all() and np.isfinite(C.lfcSE).all()

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42})
colors = {"Pupa-up": "#CC79A7", "Larva-up": "#009E73", "Non-DE": "#A7A9AC",
          "Excluded: low count": "#E69F00", "NA adjusted P": "#7E57C2", "GP-up": "#0072B2"}
WIDTH, HEIGHT = 180, 240
fig = plt.figure(figsize=(WIDTH/25.4, HEIGHT/25.4), facecolor="white")
def axmm(x, y, w, h):
    ax = fig.add_axes([x/WIDTH, y/HEIGHT, w/WIDTH, h/HEIGHT])
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.55)
    ax.set_axisbelow(True)
    return ax
def title(letter, text, x, top):
    fig.text(x/WIDTH, top/HEIGHT, letter, fontsize=11, fontweight="bold", va="top")
    fig.text((x+7)/WIDTH, top/HEIGHT, text, fontsize=9, fontweight="bold", va="top")

ax_a = axmm(15, 166, 55, 59)
title("a", "Insect inhibitor-family homologs", 4, 237)
fig.text(15/WIDTH, 229/HEIGHT, "149 candidates; 18 MEROPS families", fontsize=8)
order = list(A.category.drop_duplicates())
statuses = ["Pupa-up", "Larva-up", "Non-DE", "NA adjusted P", "Excluded: low count"]
matrix = A.pivot(index="category", columns="series", values="count").reindex(index=order, columns=statuses)
y = np.arange(len(order))
left = np.zeros(len(order))
for s in statuses:
    values = matrix[s].to_numpy()
    ax_a.barh(y, values, left=left, height=0.74, color=colors[s], edgecolor="white", linewidth=0.35)
    left += values
ax_a.set_yticks(y, order)
ax_a.set_ylim(-0.6, len(order)-0.4)
ax_a.set_xticks([0, 10, 20, 30, 40])
ax_a.set_xlabel("Candidate genes")
handles = [Patch(facecolor=colors[s], label=s) for s in statuses if matrix[s].sum()]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(6/WIDTH, 155/HEIGHT), ncol=2,
           frameon=False, fontsize=8, columnspacing=1.1, handlelength=1.1, handletextpad=0.4)

ax_b = axmm(115, 166, 57, 59)
title("b", "Larva-up inhibitor candidates", 78, 237)
fig.text(83/WIDTH, 229/HEIGHT, "Nine candidates with Pfam-family support", fontsize=8)
B = B.sort_values("legacy_display_log2FoldChange", ascending=True).reset_index(drop=True)
y = np.arange(len(B))
ax_b.hlines(y, B.plotted_log2FoldChange, 0, color="#B9B9B9", linewidth=0.6)
ax_b.scatter(B.plotted_log2FoldChange, y, s=24, color=colors["Larva-up"], edgecolor="white", linewidth=0.4, zorder=3)
ax_b.axvline(-1, color="#777777", linewidth=0.7, linestyle="--")
ax_b.set_yticks(y, B.gene_id)
ax_b.set_ylim(-0.6, len(B)-0.4)
ax_b.set_xlim(-6.75, 0.3)
ax_b.set_xticks([-6, -4, -2, 0])
ax_b.set_xlabel("log₂ fold change (pupa/larva)")
for i, row in B.iterrows():
    family = row.series.split()[0]
    label_left = row.plotted_log2FoldChange > -2.5
    ax_b.annotate(f"{row.plotted_log2FoldChange:.2f}", (row.plotted_log2FoldChange,i),
                  xytext=(-4 if label_left else 4, 0), textcoords="offset points",
                  ha="right" if label_left else "left", va="center", fontsize=8)
    ax_b.text(0.98, i, family, transform=ax_b.get_yaxis_transform(), va="center", ha="right", fontsize=8)
fig.text(82/WIDTH, 149/HEIGHT, "I01: Kazal-type; I04: serpin", fontsize=8)
fig.text(82/WIDTH, 144/HEIGHT, "I19: pacifastin", fontsize=8)

ax_c = axmm(56, 52, 105, 77)
title("c", "GP-up plant peptidase candidates", 4, 140)
fig.text(11/WIDTH, 133/HEIGHT, "Final refined-mask priority 1 + 2 classical-secretome candidates (n = 19)", fontsize=8)
C = C.sort_values("plotted_log2FoldChange", ascending=True).reset_index(drop=True)
y = np.arange(len(C))
ax_c.hlines(y, 0, C.plotted_log2FoldChange, color="#B9B9B9", linewidth=0.6)
ax_c.scatter(C.plotted_log2FoldChange, y, s=23, color=colors["GP-up"], edgecolor="white", linewidth=0.4, zorder=3)
ax_c.axvline(1, color="#777777", linewidth=0.7, linestyle="--")
ax_c.set_yticks(y, C.gene_id)
ax_c.set_ylim(-0.55, len(C)-0.45)
ax_c.set_xlim(0, 12)
ax_c.set_xticks([0, 2, 4, 6, 8, 10, 12])
ax_c.set_xlabel("log₂ fold change (GP/UP)")
for i, row in C.iterrows():
    ax_c.text(1.025, i, row.series.split(":")[0], transform=ax_c.get_yaxis_transform(),
              va="center", fontsize=8, color="#333333")
fig.text(5/WIDTH, 36/HEIGHT, "MEROPS: A01 aspartic; C01 papain-like; C13 legumain/VPE; M10 metallo;", fontsize=8)
fig.text(5/WIDTH, 31/HEIGHT, "S08 subtilisin-like; S10 serine carboxypeptidase. All points: GP-up.", fontsize=8)

ax_d = axmm(79, 8, 88, 13)
title("d", "Insect inhibitor-set stage bias", 4, 28)
D = D.reset_index(drop=True)
y = np.arange(len(D))
ax_d.errorbar(D.odds_ratio, y,
    xerr=[D.odds_ratio-D.confidence_low, D.confidence_high-D.odds_ratio], fmt="o", markersize=4,
    color=colors["Pupa-up"], ecolor="#555555", elinewidth=0.8, capsize=3)
ax_d.axvline(1, linestyle="--", linewidth=0.7, color="#777777")
ax_d.set_xscale("log", base=2)
ax_d.set_xlim(0.8, 16)
ax_d.set_xticks([1, 2, 4, 8, 16], ["1", "2", "4", "8", "16"])
ax_d.set_yticks(y, ["All curated A+B", "Tier A named-reference"])
ax_d.set_ylim(-0.5, 1.8)
for i, row in D.iterrows():
    ax_d.text(15.4, i+0.3, f"OR {row.odds_ratio:.2f}; FDR {row.padj:.2g}",
              fontsize=8, va="bottom", ha="right")
fig.text(79/WIDTH, 1/HEIGHT, "Odds ratio (95% CI); >1: pupa-up enrichment", fontsize=8)

stem = args.outdir / "Figure_10_cross_system_MEROPS"
fig.savefig(stem.with_suffix(".pdf"), metadata={"Creator": "Figure_10_plotting_script.py"})
fig.savefig(stem.with_suffix(".png"), dpi=300)
plt.close(fig)
print("Figure 10: source rows 120; panels A/B/C/D 90/9/19/2; insect homologs 149; inhibitor families 18; larva shortlist 9; GP-up priority candidates 19.")
