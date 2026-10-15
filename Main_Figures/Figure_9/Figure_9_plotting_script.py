#!/usr/bin/env python3
"""Reproduce Figure 9 from the frozen refined-mask plant-defense overlay.

Revision 2026-10-09: retain the ten candidates shown in the manuscript
and supplied PNG. Correct the prior script/PDF substitution of DN3635
for DN27036. This script only summarizes existing statuses and plots
frozen values; it does not fit models, annotate, or select new candidates.
Requires Python 3, matplotlib, numpy, and pandas.
Run from any directory: python /path/to/Figure_9_plotting_script.py
"""

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "plant_refined_mask_defense_candidate_gene_summary.tsv"

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--outdir", type=Path, default=SCRIPT_DIR)
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)

data = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
assert len(data) == 806
assert data["gene_id"].nunique() == 806
assert data["candidate_membership_changed"].str.upper().eq("FALSE").all()

status_order = [
    "GP_up",
    "UP_up",
    "non_DE",
    "NA_padj_not_testable",
    "excluded_by_refined_sequence_mask",
    "not_in_refined_mask_host_focused_tested_universe",
]
status_labels = {
    "GP_up": "GP-up",
    "UP_up": "UP-up",
    "non_DE": "Non-DE",
    "NA_padj_not_testable": "NA adjusted P",
    "excluded_by_refined_sequence_mask": "Excluded by refined mask",
    "not_in_refined_mask_host_focused_tested_universe": "Outside final tested universe",
}
status_colors = {
    "GP_up": "#0072B2",
    "UP_up": "#D55E00",
    "non_DE": "#A7A9AC",
    "NA_padj_not_testable": "#7E57C2",
    "excluded_by_refined_sequence_mask": "#E69F00",
    "not_in_refined_mask_host_focused_tested_universe": "#D7DBE0",
}

status_counts = (
    data["refined_mask_host_status"].value_counts().reindex(status_order, fill_value=0).astype(int)
)
assert status_counts.to_dict() == {
    "GP_up": 86,
    "UP_up": 28,
    "non_DE": 518,
    "NA_padj_not_testable": 14,
    "excluded_by_refined_sequence_mask": 3,
    "not_in_refined_mask_host_focused_tested_universe": 157,
}

group_order = [
    "PR_protein",
    "jasmonate",
    "salicylate",
    "peptidase",
    "peptidase_inhibitor",
]
group_labels = {
    "PR_protein": "PR proteins",
    "jasmonate": "Jasmonate-associated",
    "salicylate": "Salicylate-associated",
    "peptidase": "Peptidases",
    "peptidase_inhibitor": "Peptidase inhibitors",
}

memberships = data["defense_groups"].str.split(";").apply(
    lambda values: {value.strip() for value in values}
)
group_de_counts = pd.DataFrame(index=group_order, columns=["GP_up", "UP_up"], data=0)
for group in group_order:
    member_mask = memberships.apply(lambda groups: group in groups)
    counts = data.loc[member_mask, "refined_mask_host_status"].value_counts()
    group_de_counts.loc[group, "GP_up"] = int(counts.get("GP_up", 0))
    group_de_counts.loc[group, "UP_up"] = int(counts.get("UP_up", 0))

expected_group_de = pd.DataFrame(
    [[39, 8], [5, 7], [2, 0], [40, 13], [8, 0]],
    index=group_order,
    columns=["GP_up", "UP_up"],
)
assert group_de_counts.equals(expected_group_de)

selected = {
    "TRINITY_DN99535_c0_g1": "PR6/Kunitz inhibitor",
    "TRINITY_DN8134_c0_g1": "PR6/Kunitz inhibitor",
    "TRINITY_DN60637_c0_g1": "PR6/Kunitz inhibitor",
    "TRINITY_DN4958_c0_g1": "PR6/Kunitz inhibitor",
    "TRINITY_DN11853_c0_g1": "PR9 peroxidase",
    "TRINITY_DN8699_c0_g1": "Peptidase",
    "TRINITY_DN11522_c0_g1": "Peptidase",
    "TRINITY_DN27036_c0_g1": "Peptidase",
    "TRINITY_DN13180_c0_g1": "Peptidase",
    "TRINITY_DN520_c1_g1": "Peptidase",
}
selected_df = data.loc[data["gene_id"].isin(selected)].copy()
assert len(selected_df) == len(selected)
assert selected_df["refined_mask_host_status"].eq("GP_up").all()
assert selected_df["signalp_prediction"].eq("SP").all()
assert selected_df["has_TMR"].str.upper().eq("FALSE").all()
selected_df["candidate_class"] = selected_df["gene_id"].map(selected)
for column in [
    "refined_mask_host_baseMean",
    "refined_mask_host_log2FoldChange",
    "refined_mask_host_padj",
]:
    selected_df[column] = pd.to_numeric(selected_df[column], errors="raise")
selected_df = selected_df.sort_values("refined_mask_host_log2FoldChange", ascending=True)

# Save the exact data displayed in the figure.
source_rows = []
for status in status_order:
    source_rows.append(
        {
            "panel": "A",
            "category": status_labels[status],
            "series": "Unique defense candidates",
            "count": int(status_counts[status]),
            "gene_id": "",
            "candidate_class": "",
            "baseMean": np.nan,
            "log2FoldChange": np.nan,
            "padj": np.nan,
        }
    )
for group in group_order:
    for status in ["GP_up", "UP_up"]:
        source_rows.append(
            {
                "panel": "B",
                "category": group_labels[group],
                "series": status_labels[status],
                "count": int(group_de_counts.loc[group, status]),
                "gene_id": "",
                "candidate_class": "",
                "baseMean": np.nan,
                "log2FoldChange": np.nan,
                "padj": np.nan,
            }
        )
for _, row in selected_df.iterrows():
    source_rows.append(
        {
            "panel": "C",
            "category": "Selected GP-up candidates",
            "series": row["candidate_class"],
            "count": np.nan,
            "gene_id": row["gene_id"],
            "candidate_class": row["candidate_class"],
            "baseMean": data.set_index("gene_id").at[row["gene_id"], "refined_mask_host_baseMean"],
            "log2FoldChange": data.set_index("gene_id").at[row["gene_id"], "refined_mask_host_log2FoldChange"],
            "padj": data.set_index("gene_id").at[row["gene_id"], "refined_mask_host_padj"],
            "expression_status": row["refined_mask_host_status"],
            "signalp_prediction": row["signalp_prediction"],
            "has_TMR": row["has_TMR"],
            "MEROPS_families": row["MEROPS_families"],
        }
    )
source_df = pd.DataFrame(source_rows)
source_df["count"] = source_df["count"].astype("Int64")
source_df.to_csv(
    args.outdir / "Figure_9_source_data.tsv", sep="\t", index=False
)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 9.5,
        "axes.titlesize": 11.5,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

fig = plt.figure(figsize=(14.2, 10.0), constrained_layout=False)
grid = fig.add_gridspec(
    2, 2, height_ratios=[0.84, 1.16], width_ratios=[1.0, 1.0], hspace=0.53, wspace=0.39
)
ax_a = fig.add_subplot(grid[0, 0])
ax_b = fig.add_subplot(grid[0, 1])
ax_c = fig.add_subplot(grid[1, :])

# Panel A: mutually exclusive final status across all 806 unique candidates.
y = np.arange(len(status_order))
bars = ax_a.barh(
    y,
    status_counts.values,
    color=[status_colors[status] for status in status_order],
    edgecolor="#3F3F3F",
    linewidth=0.55,
)
ax_a.set_yticks(y, [status_labels[status] for status in status_order])
ax_a.invert_yaxis()
ax_a.set_xlabel("Unique candidate genes")
ax_a.set_title("Final expression status of the frozen inventory", loc="left", fontweight="bold")
ax_a.grid(axis="x", color="#E3E6E8", linewidth=0.7)
ax_a.set_axisbelow(True)
ax_a.set_xlim(0, 580)
for bar, value in zip(bars, status_counts.values):
    ax_a.text(
        bar.get_width() + 8,
        bar.get_y() + bar.get_height() / 2,
        f"{value:,}",
        va="center",
        fontsize=9,
    )

# Panel B: DE direction within each overlapping defense group.
y = np.arange(len(group_order))
height = 0.35
gp = ax_b.barh(
    y - height / 2,
    group_de_counts["GP_up"],
    height=height,
    color=status_colors["GP_up"],
    edgecolor="#3F3F3F",
    linewidth=0.5,
    label="GP-up",
)
up = ax_b.barh(
    y + height / 2,
    group_de_counts["UP_up"],
    height=height,
    color=status_colors["UP_up"],
    edgecolor="#3F3F3F",
    linewidth=0.5,
    label="UP-up",
)
ax_b.set_yticks(y, [group_labels[group] for group in group_order])
ax_b.invert_yaxis()
ax_b.set_xlabel("Differentially expressed candidate genes")
ax_b.set_title("Direction of differential expression by defense group", loc="left", fontweight="bold")
ax_b.grid(axis="x", color="#E3E6E8", linewidth=0.7)
ax_b.set_axisbelow(True)
ax_b.set_xlim(0, 46)
ax_b.legend(loc="lower right", frameon=False)
for collection in [gp, up]:
    for bar in collection:
        value = int(bar.get_width())
        if value > 0:
            ax_b.text(value + 0.7, bar.get_y() + bar.get_height() / 2, str(value), va="center", fontsize=8.5)
ax_b.text(
    0.0,
    -0.20,
    "Groups overlap; counts must not be summed across categories.",
    transform=ax_b.transAxes,
    fontsize=8.3,
    color="#555555",
)

# Panel C: selected GP-up candidates reported in Section 3.8.
class_colors = {
    "PR6/Kunitz inhibitor": "#009E73",
    "PR9 peroxidase": "#7E57C2",
    "Peptidase": "#56B4E9",
}
y = np.arange(len(selected_df))
for candidate_class in class_colors:
    mask = selected_df["candidate_class"].eq(candidate_class).to_numpy()
    ax_c.scatter(
        selected_df.loc[mask, "refined_mask_host_log2FoldChange"],
        y[mask],
        s=75,
        color=class_colors[candidate_class],
        edgecolor="white",
        linewidth=0.8,
        label=candidate_class,
        zorder=3,
    )
ax_c.axvline(1, linestyle="--", linewidth=1.0, color="#777777")
ax_c.set_yticks(y, selected_df["gene_id"])
ax_c.set_xlabel("log₂ fold change (GP / UP)")
ax_c.set_title(
    "Selected GP-up candidates with classical secretory-pathway evidence",
    loc="left",
    fontweight="bold",
)
ax_c.grid(axis="x", color="#E3E6E8", linewidth=0.7)
ax_c.set_axisbelow(True)
ax_c.set_xlim(0, 14.1)
ax_c.legend(loc="lower right", frameon=False, ncol=3)
for index, row in selected_df.reset_index(drop=True).iterrows():
    ax_c.text(
        row["refined_mask_host_log2FoldChange"] + 0.18,
        index,
        f'{row["refined_mask_host_log2FoldChange"]:.2f}',
        va="center",
        fontsize=8.3,
        color="#333333",
    )
ax_c.text(
    0.99,
    -0.15,
    "All displayed candidates were SignalP-positive and lacked a predicted transmembrane region.",
    transform=ax_c.transAxes,
    ha="right",
    fontsize=8.3,
    color="#555555",
)

for label, ax in zip("ABC", [ax_a, ax_b, ax_c]):
    ax.text(-0.14 if label != "C" else -0.065, 1.08, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

fig.suptitle(
    "Focused plant-defense candidates in galled and ungalled petioles",
    x=0.07,
    ha="left",
    y=0.985,
    fontsize=15,
    fontweight="bold",
)
fig.subplots_adjust(left=0.20, right=0.975, top=0.92, bottom=0.085)

fig.savefig(args.outdir / "Figure_9_plant_defense_candidates.pdf", bbox_inches="tight")
fig.savefig(
    args.outdir / "Figure_9_plant_defense_candidates.png",
    dpi=300,
    bbox_inches="tight",
    pil_kwargs={"compress_level": 6},
)
plt.close(fig)

print("Figure 9 created from the frozen 806-gene plant-defense overlay")
