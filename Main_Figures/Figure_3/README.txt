Figure 3: insect volcano. All 13,919 source rows are included; 13,915 rows with available adjusted P values are plotted.
Four unavailable-P rows are retained in source data with plotted = FALSE. Non-DE denotes genes with available adjusted P values not meeting both DEG cutoffs.
The y axis is capped at 50, with 237 points exceeding 50 displayed at the cap; actual values remain in source data.
Positive x values denote Pupa-up and negative x values denote Larva-up.
The source table preserves gene IDs, baseMean, unshrunk log2FoldChange, lfcSE, Wald statistic, raw P, and adjusted P directly from the validated result.
Original numeric strings are retained; explicit Python float conversion avoids pandas fast-parser loss of very small fixed-point P values.

Rebuild: python plot.py
Dependencies: Python 3, matplotlib, numpy.
The script resolves all inputs relative to itself and may be run from any working directory.
PDF width: 180 mm; PNG resolution: 300 dpi. Smallest text: 10 pt.
All edits are display-only; no differential-expression analysis or PCA was rerun.
