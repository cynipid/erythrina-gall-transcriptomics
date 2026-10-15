Figure 10 reproduction and revision notes

Run: python Figure_10_plotting_script.py
The script resolves the adjacent Figure_10_source_data.tsv relative to its own
location and can be run from any working directory. Requires Python 3, pandas,
numpy, and matplotlib. --source and --outdir may override those locations.

Output: fixed 180 x 240 mm vector PDF and 300 dpi PNG. Minimum text size is 8 pt.
No MEROPS search, differential-expression test, enrichment test, candidate
selection, or shrinkage is rerun. The same 120 source rows are retained:
panel a 90 rows (18 families x 5 expression categories), panel b 9 genes,
panel c 19 genes, and panel d 2 candidate-set tests.

Panel a: 149 inhibitor-family homologs in the curated A+B set, stacked by status.
The display label 'Not a primary DEG' is standardized to 'Non-DE'.
Panel b: the same 9 larva-up inhibitors with expected Pfam-family support.
The horizontal axis now shows the original unshrunk log2 fold change (pupa/larva),
so larva-up values are negative and the dashed threshold is -1. Full source
precision is retained. Family codes identify I01 (Kazal-type), I04 (serpin), and
I19 (pacifastin); all points are green to encode larva-up consistently.
Panel c: the same 19 final refined-mask GP-up priority 1 + 2 peptidase candidates.
The horizontal axis retains unshrunk log2 fold change (GP/UP), threshold +1.
All points are blue for GP-up; MEROPS group codes are shown to the right.
Panel d: the same odds ratios, confidence intervals, and adjusted P values from
the existing Fisher exact tests. Pink points indicate pupa-up enrichment.

Color semantics: larva-up #009E73; pupa-up #CC79A7; GP-up #0072B2;
Non-DE #A7A9AC. Non-DE means failure of the combined primary criterion (adjusted
P < 0.01 AND absolute log2 fold change >= 1), not absence of expression.

Source columns explicitly preserve sign provenance:
- legacy_display_log2FoldChange: value in the previous figure's source data.
- original_log2FoldChange / original_contrast: the original DESeq2 orientation.
- plotted_log2FoldChange / plotted_contrast: values and orientation now plotted.
- log2FoldChange: matches original_log2FoldChange for direct conventional use.
- coefficient_type: unshrunk for panels b and c.
- lfcSE: frozen DESeq2 standard error for all 9 insect and 19 plant candidates.
For panel b only, original = plotted = -legacy_display. Panel c is unchanged.
All other panels have no fold-change value. No coefficient was truncated.

Source verification: original insect coefficients and padj match the curated
MEROPS candidate table and frozen insect DESeq2 all-gene table. Plant coefficients
and padj match the final refined-mask defense overlay and all-gene table.
The lfcSE range is 0.201656462238479 to 0.938653390267582 for insect candidates,
and 0.314884107286783 to 1.80548639633566 for plant candidates. These standard
errors are supplied in the source data; the plotted horizontal stems in b and c
connect coefficients to zero and are not confidence intervals. Panel d displays
the original 95% confidence intervals for odds ratios.
