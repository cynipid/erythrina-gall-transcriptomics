Figure 6: display revision, 26 October 2026

The 180 x 240 mm canvas retains four panels. Panel c shows the same nineteen full
Trinity identifiers at 10 pt, with about 4.45 mm row spacing. Panel-b identifiers
are 9 pt. Spacing was adjusted; source values, signs, order, colors, axis ranges,
odds ratios and confidence intervals are unchanged.
Figure_6_source_data.tsv retains all 120 frozen rows (A/B/C/D: 90/9/19/2).

Rebuild: python plot.py
Write elsewhere: python plot.py --outdir /path/to/test_output
The script reads only its adjacent source TSV; no expression, enrichment or
MEROPS analysis is rerun. PDF and 300-dpi PNG are generated together.
