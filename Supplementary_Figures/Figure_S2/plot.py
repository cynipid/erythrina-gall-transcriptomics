#!/usr/bin/env python3
"""Plot only the supplied full-precision correlation matrix.

Run: python plot.py
Dependencies: Python 3, numpy, pandas, matplotlib.
The final production PDF/PNG use this portable Matplotlib rendering.
Every coefficient and sample label is preserved from the accepted matrix.
The R-generated artwork was cosmetically refreshed for script/output consistency.
No model, transformation, PCA or correlation calculation is performed here.
"""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

import argparse

HERE = Path(__file__).resolve().parent
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--outdir', type=Path, default=None, help='Output directory; source inputs remain beside this script.')
_args = _parser.parse_args()
OUTDIR = _args.outdir.expanduser().resolve() if _args.outdir is not None else HERE
OUTDIR.mkdir(parents=True, exist_ok=True)
CONFIG = {'source': 'Supplementary_Figure_S2_source_data.tsv', 'stem': 'Supplementary_Figure_S2', 'labels': ['L1', 'L2', 'L3', 'P1', 'P2', 'P3'], 'colors': ['#2166AC', 'white', '#B2182B'], 'limits': [0.6, 1], 'ticks': [0.6, 0.7, 0.8, 0.9, 1], 'rotation': 0}
source = HERE / CONFIG['source']
matrix = pd.read_csv(source, sep='\t', index_col=0)
labels = CONFIG['labels']
assert matrix.shape == (6, 6)
assert matrix.index.tolist() == labels and matrix.columns.tolist() == labels
values = matrix.to_numpy(dtype=float)
assert np.isfinite(values).all() and np.allclose(values, values.T, atol=1e-14)
assert np.allclose(np.diag(values), 1, atol=1e-14)
limits = CONFIG['limits'] or [float(values.min()), float(values.max())]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
fig = plt.figure(figsize=(6.6, 5.8))
ax = fig.add_axes([0.075, 0.13, 0.735, 0.81])
cmap = LinearSegmentedColormap.from_list('retained_palette', CONFIG['colors'], N=1024)
mesh = ax.pcolormesh(np.arange(7), np.arange(7), values, cmap=cmap,
                     vmin=limits[0], vmax=limits[1], edgecolors='white', linewidth=1.3)
ax.set_aspect('equal'); ax.invert_yaxis()
ax.set_xticks(np.arange(6) + .5, labels=labels,
              rotation=CONFIG['rotation'], ha='right' if CONFIG['rotation'] else 'center')
ax.set_yticks(np.arange(6) + .5, labels=labels)
ax.tick_params(length=0, pad=7)
for edge in ax.spines.values(): edge.set_visible(False)
for i in range(6):
    for j in range(6):
        ax.text(j + .5, i + .5, f'{values[i,j]:.2f}', ha='center', va='center', fontsize=12)
cax = fig.add_axes([0.86, 0.375, 0.025, 0.255])
bar = fig.colorbar(mesh, cax=cax, ticks=CONFIG['ticks'])
bar.outline.set_visible(False); bar.ax.tick_params(length=0, pad=5, labelsize=10)
bar.ax.set_title('Pearson r', fontsize=11, loc='center', pad=13)
for ext in ['pdf','png']:
    fig.savefig(OUTDIR / (CONFIG['stem'] + '.' + ext), dpi=400, facecolor='white',
                metadata={'CreationDate': None, 'ModDate': None} if ext == 'pdf' else None)
plt.close(fig)
audit = {'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
         'sample_order': labels, 'matrix_shape': [6,6],
         'minimum': float(values.min()), 'maximum': float(values.max()),
         'all_cells_displayed_to_two_decimal_places': True,
         'numeric_matrix_modified': False, 'analysis_rerun': False,
         'rendering_note': 'Final production artwork uses this portable Matplotlib script; source matrix and sample labels unchanged.'}
(OUTDIR / 'render_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps(audit, indent=2))
