#!/usr/bin/env python3
"""Reproduce Figure 5 using only the adjacent frozen source TSV.

Requires Python 3, pandas, numpy, and matplotlib. Run from any working directory:
    python /path/to/Figure_5_plotting_script.py
No differential-expression analysis or candidate selection is rerun.
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
parser.add_argument("--source", type=Path, default=HERE / "Figure_5_source_data.tsv")
parser.add_argument("--outdir", type=Path, default=HERE)
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)
source = pd.read_csv(args.source, sep="\t")
assert source.groupby("panel").size().to_dict() == {"A": 5, "B": 5, "C": 25, "D": 14}
A, B, C, D = [source[source.panel.eq(p)].copy() for p in "ABCD"]
assert int(A.value.sum()) == 13455 and int(B.value.sum()) == 2039
assert int(C.value.sum()) == 66 and D.gene_id.nunique() == 14
assert D.log2FoldChange.le(-1).all() and D.padj.lt(0.01).all()
assert np.isfinite(D.lfcSE).all() and D.lfcSE.gt(0).all()

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42})
colors = {"Pupa-up": "#CC79A7", "Larva-up": "#009E73", "Non-DE": "#A7A9AC",
          "Excluded: low count": "#E69F00", "NA adjusted P": "#7E57C2"}
family_short = {"GH12 cellulase/endoglucanase": "GH12 cellulase/endoglucanase",
    "Pectin/pectate lyase": "Pectin/pectate lyase",
    "Canonical GH18 chitinase": "Canonical GH18 chitinase",
    "GH19 class-I chitinase": "GH19 class-I chitinase", "CAP-domain protein": "CAP-domain protein"}
family_tags = {"GH12 cellulase/endoglucanase": "GH12", "Pectin/pectate lyase": "Lyase",
    "Canonical GH18 chitinase": "GH18", "GH19 class-I chitinase": "GH19", "CAP-domain protein": "CAP"}
WIDTH, HEIGHT = 180, 238
fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor="white")
def axmm(x, y, w, h):
    ax = fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.55)
    ax.set_axisbelow(True)
    return ax
def title(letter, text, top):
    fig.text(4 / WIDTH, top / HEIGHT, letter, fontsize=11, fontweight="bold", va="top")
    fig.text(11 / WIDTH, top / HEIGHT, text, fontsize=9.5, fontweight="bold", va="top")

ax_a = axmm(65, 198, 107, 29)
title("a", "SignalP/DeepTMHMM classification", 234)
y = np.arange(len(A))
bars = ax_a.barh(y, A.value, color=["#B8BDC5", "#6C7A89", "#444444", "#E69F00", "#7E57C2"], height=0.67)
ax_a.set_yticks(y, A.category)
ax_a.invert_yaxis()
ax_a.set_xlim(0, 10300)
ax_a.set_xticks([0, 2500, 5000, 7500, 10000], ["0", "2,500", "5,000", "7,500", "10,000"])
ax_a.set_xlabel("Representative proteins")
for bar, count in zip(bars, A.value):
    ax_a.text(count + 100, bar.get_y() + bar.get_height()/2, f"{int(count):,}", va="center", fontsize=8)

ax_b = axmm(65, 152, 107, 29)
title("b", "Expression status of classical-secretome candidates", 188)
y = np.arange(len(B))
bars = ax_b.barh(y, B.value, color=[colors[x] for x in B.category], height=0.67)
ax_b.set_yticks(y, B.category)
ax_b.invert_yaxis()
ax_b.set_xlim(0, 870)
ax_b.set_xticks([0, 200, 400, 600, 800])
ax_b.set_xlabel("Classical-secretome candidates")
for bar, count in zip(bars, B.value):
    ax_b.text(count + 8, bar.get_y() + bar.get_height()/2, f"{int(count):,}", va="center", fontsize=8)

ax_c = axmm(65, 106, 107, 29)
title("c", "Expression status within focused families", 142)
family_order = list(C.category.drop_duplicates())
status_order = list(B.category)
matrix = C.pivot(index="category", columns="series", values="value").reindex(index=family_order, columns=status_order)
y = np.arange(len(matrix))
left = np.zeros(len(matrix))
for status in status_order:
    values = matrix[status].to_numpy()
    ax_c.barh(y, values, left=left, height=0.73, color=colors[status], edgecolor="white", linewidth=0.5)
    for j, (n, x) in enumerate(zip(values, left)):
        if n > 0:
            ax_c.text(x+n/2, j, str(int(n)), ha="center", va="center", color="#111111", fontsize=8)
    left += values
ax_c.set_yticks(y, [family_short[x] for x in family_order])
ax_c.invert_yaxis()
ax_c.set_xlim(0, 17)
ax_c.set_xticks([0, 4, 8, 12, 16])
ax_c.set_xlabel("Candidate genes")
handles = [Patch(facecolor=colors[s], label=s) for s in status_order if matrix[s].sum()]
fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.53, 91/HEIGHT), ncol=4,
           frameon=False, fontsize=8, handlelength=1.1, columnspacing=1.2, handletextpad=0.4)

ax_d = axmm(65, 17, 107, 62)
title("d", "Larva-up focused-family classical-secretome candidates", 85)
D = D.sort_values("baseMean", ascending=True).reset_index(drop=True)
y = np.arange(len(D))
ax_d.scatter(D.baseMean, y, s=22, color=colors["Larva-up"], linewidths=0.45, edgecolors="white", zorder=3)
ax_d.set_xscale("log")
ax_d.set_xticks([100, 1000, 10000], ["100", "1,000", "10,000"])
ax_d.set_yticks(y, [f"{r.gene_id}  [{family_tags[r.category]}]" for _, r in D.iterrows()])
ax_d.set_ylim(-0.55, len(D)-0.45)
ax_d.set_xlim(D.baseMean.min()/1.3, D.baseMean.max()*2.2)
ax_d.set_xlabel("Mean normalized count (baseMean; log scale)")
for i, row in D.iterrows():
    ax_d.annotate(f"{row.log2FoldChange:.2f}", (row.baseMean, i), xytext=(4,0),
                  textcoords="offset points", va="center", fontsize=8)
fig.text(65/WIDTH, 2/HEIGHT, "Numbers: unshrunk log₂ fold change (pupa/larva).", fontsize=8)

stem = args.outdir / "Figure_5_insect_secretome_candidate_families"
fig.savefig(stem.with_suffix(".pdf"), metadata={"Creator": "Figure_5_plotting_script.py"})
fig.savefig(stem.with_suffix(".png"), dpi=300)
plt.close(fig)
print("Figure 5: source rows 49; panels A/B/C/D 5/5/25/14; proteins 13,455; classical-secretome candidates 2,039; focused genes 66; larva-up shortlist 14.")
