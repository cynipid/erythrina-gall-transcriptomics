Figure 5 reproduction and revision notes

Run: python Figure_5_plotting_script.py
The script resolves the adjacent Figure_5_source_data.tsv relative to its own
location and can be run from any working directory. Requires Python 3, pandas,
numpy, and matplotlib. --source and --outdir may override those locations.

Output: fixed 180 x 238 mm vector PDF and 300 dpi PNG. Minimum text size is 8 pt.
The four panels are retained, now stacked to preserve readable family and gene
labels. No enrichment, differential-expression test, or candidate selection is
rerun. All 49 original plotted source rows and their values are retained.

Panel a: 13,455 representative proteins in five prediction classes.
Panel b: 2,039 classical-secretome candidates: 597 pupa-up, 544 larva-up,
790 Non-DE, 105 excluded for low counts, and 3 with NA adjusted P. The horizontal
axis now includes 790 fully; the previous plot's 660 limit clipped this bar.
Panel c: 66 focused-family genes; all 25 family-by-status rows are preserved,
including zero values. Counts of one are now labeled as well.
Panel d: the same 14 larva-up classical-secretome candidates. Uniform green
points indicate larva-up; family abbreviations appear beside the unchanged full
gene IDs. Point sizes no longer encode effect size. Numeric labels retain the
unshrunk log2 fold change in the original pupa/larva orientation.

Color semantics: larva-up #009E73; pupa-up #CC79A7; Non-DE #A7A9AC.
Prediction classes in panel a have their own non-stage colors.
Non-DE means the gene fails the combined primary DEG criterion (adjusted P <
0.01 AND absolute log2 fold change >= 1), not absence of expression.

Uncertainty: the source TSV adds lfcSE for all 14 panel-d candidates, copied from
the frozen insect all-gene DESeq2 results. No shrinkage was performed. lfcSE is
the standard error of the unshrunk coefficient; it is not displayed as a count
error bar because the horizontal axis in panel d is baseMean. baseMean and padj
remain available. lfcSE range: 0.335365058965921 to 1.83448773859676.
