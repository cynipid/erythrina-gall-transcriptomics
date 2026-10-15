Figure 7: plant volcano using the final refined mask. All 24,734 rows are plotted: 2,147 GP-up, 1,383 UP-up, 21,032 Non-DE, and 172 Not tested.
Not tested denotes unavailable adjusted P values; all 172 also lack raw P values. These observations are shown as purple crosses at y = 0 only for visualization, with no P-value imputation.
The y axis is not capped; minimum adjusted P = 2.3355264667008e-37.
Extreme unshrunk log2 fold changes can be inspected alongside baseMean and lfcSE in the source table; no shrinkage sensitivity analysis is claimed.
Positive x values denote GP-up and negative x values denote UP-up.
The source table preserves gene IDs, baseMean, unshrunk log2FoldChange, lfcSE, Wald statistic, raw P, and adjusted P directly from the validated result.
Original numeric strings are retained; explicit Python float conversion avoids pandas fast-parser loss of very small fixed-point P values.

Rebuild: python plot.py
Dependencies: Python 3, matplotlib, numpy.
The script resolves all inputs relative to itself and may be run from any working directory.
PDF width: 180 mm; PNG resolution: 300 dpi. Smallest text: 10 pt.
All edits are display-only; no differential-expression analysis or PCA was rerun.
