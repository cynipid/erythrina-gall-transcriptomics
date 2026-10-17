#!/usr/bin/env python3
"""Render an MA plot from frozen DESeq2 result rows; no model is fitted.

Requirements: Python 3, numpy, pandas and matplotlib.
Run: python plot.py
All inputs and outputs resolve relative to this script's directory.
The source TSV preserves the archived volcano-source bytes. Volcano-specific
display fields are ignored: MA coordinates use original baseMean and unshrunk
log2FoldChange. Classification is read verbatim from status, not recalculated.
"""
from pathlib import Path
from collections import Counter
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import argparse

HERE = Path(__file__).resolve().parent
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--outdir', type=Path, default=None, help='Output directory; source inputs remain beside this script.')
_args = _parser.parse_args()
OUTDIR = _args.outdir.expanduser().resolve() if _args.outdir is not None else HERE
OUTDIR.mkdir(parents=True, exist_ok=True)
CONFIG = {'source': 'Supplementary_Figure_S8_source_data.tsv', 'stem': 'Supplementary_Figure_S8', 'direction': 'GP / UP', 'negative': 'UP-up', 'positive': 'GP-up', 'counts': {'GP-up': 2147, 'UP-up': 1383, 'Non-DE': 21032, 'Not tested': 172}, 'colors': {'GP-up': '#0072B2', 'UP-up': '#D55E00', 'Non-DE': '#BDBDBD'}}
source = HERE / CONFIG['source']
data = pd.read_csv(source, sep='\t')
expected = CONFIG['counts']
observed = dict(Counter(data['status']))
assert observed == expected, (observed, expected)
assert data['gene_id'].is_unique, 'Duplicate gene identifiers'
assert np.isfinite(data['baseMean']).all() and (data['baseMean'] >= 0).all()
assert np.isfinite(data['log2FoldChange']).all(), 'No fabricated coordinates permitted'
assert data.loc[data['status'].eq('Not tested'), 'padj'].isna().all()
assert data.loc[~data['status'].eq('Not tested'), 'padj'].notna().all()

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.labelsize': 11, 'axes.linewidth': 0.8,
    'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
fig, ax = plt.subplots(figsize=(7.2, 5.4))
fig.subplots_adjust(left=0.12, right=0.985, bottom=0.135, top=0.82)

# Draw complete frozen rows. Not-tested genes retain their actual finite LFCs.
draw_order = ['Non-DE', CONFIG['negative'], CONFIG['positive'], 'Not tested']
for status in draw_order:
    subset = data.loc[data['status'].eq(status)]
    x = np.log10(subset['baseMean'].to_numpy() + 1.0)
    y = subset['log2FoldChange'].to_numpy()
    if status == 'Not tested':
        ax.scatter(x, y, marker='x', s=19, color='#666666',
                   linewidths=0.85, alpha=0.9, zorder=4)
    else:
        ax.scatter(x, y, marker='o', s=7,
                   color=CONFIG['colors'][status], edgecolors='none',
                   alpha=0.45 if status == 'Non-DE' else 0.65,
                   zorder=1 if status == 'Non-DE' else 2)

ax.axhline(0, color='#333333', linewidth=0.75, zorder=0)
for threshold in [-1, 1]:
    ax.axhline(threshold, color='#777777', linewidth=0.7,
               linestyle=(0, (5, 4)), zorder=0)
ax.set_xlabel(r'$\log_{10}(\mathrm{baseMean}+1)$')
ax.set_ylabel(r'Unshrunk $\log_2$ fold change (' + CONFIG['direction'] + ')')
ax.spines[['top', 'right']].set_visible(False)
ax.tick_params(axis='both', direction='out', length=3.5, width=0.75)
ax.margins(x=0.025, y=0.065)

legend_order = [CONFIG['negative'], CONFIG['positive'], 'Non-DE', 'Not tested']
handles = []
for status in legend_order:
    if status == 'Not tested':
        handle = Line2D([], [], marker='x', color='#666666', linestyle='None',
                        markersize=5, markeredgewidth=0.9)
    else:
        handle = Line2D([], [], marker='o', color=CONFIG['colors'][status],
                        linestyle='None', markersize=4.5, markeredgewidth=0)
    handle.set_label(f'{status} (n = {expected[status]:,})')
    handles.append(handle)
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.54, 0.99),
           ncol=2, frameon=False, handletextpad=0.5, columnspacing=2.1,
           labelspacing=0.55, fontsize=9.5)

for extension in ['pdf', 'png']:
    fig.savefig(OUTDIR / (CONFIG['stem'] + '.' + extension),
                dpi=300, facecolor='white')
plt.close(fig)

audit = {
    'source': source.name,
    'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
    'rows': len(data), 'status_counts': observed,
    'x_coordinate': 'log10(baseMean + 1)',
    'y_coordinate': 'unshrunk log2FoldChange',
    'positive_direction': CONFIG['positive'],
    'negative_direction': CONFIG['negative'],
    'not_tested_coordinates': 'Actual finite log2FoldChange; no artificial zeros',
    'classification': 'Read unchanged from archived status column',
    'statistical_analysis_rerun': False,
    'output_png_dpi': 300,
}
(OUTDIR / 'render_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps(audit, indent=2))
