Supplementary Figure S8 - MA plot

Run python plot.py. Requirements: Python 3, numpy, pandas and matplotlib.
The frozen source TSV is byte-identical to the authoritative archived volcano
source table; volcano-specific display columns are ignored.
The x-axis uses log10(baseMean + 1), and the y-axis uses original unshrunk
log2FoldChange. All statuses are read unchanged from the source.
Not-tested genes have unavailable adjusted P values and are plotted as crosses
at actual finite fold-change coordinates. No artificial fold-change values
are introduced. No DESeq2 fitting, filtering, shrinkage, or analysis rerun occurs.
The primary DEG rule remains adjusted P < 0.01 and absolute log2FC >= 1.
