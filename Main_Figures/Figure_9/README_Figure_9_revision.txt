FIGURE 9 CORRECTION 2026-10-09

Purpose
Reconcile the Figure 9 PDF and plotting script with the ten candidates already
shown in the submitted manuscript and supplied PNG. This package addresses the
Figure 9 identifier mismatch and missing source data from the revision list.

Decision and evidence
Retain TRINITY_DN27036_c0_g1 in panel C. The manuscript Results paragraph beginning
"The GP-up peptidase subset included candidates from several catalytic families"
and its embedded Figure 9 both identify this gene. The frozen input table gives:
  final status: GP_up
  log2FoldChange: 10.9964531687837 (displayed as 11.00)
  adjusted P value: 0.0000000586696732521244
  SignalP: SP
  transmembrane region: FALSE
  MEROPS family: C01A
The received PDF and script instead included TRINITY_DN3635_c0_g1, a different
valid GP-up candidate (log2FoldChange 10.3145261629095). That alternative has not
been removed from the 806-gene inventory. Only the figure selection was aligned
with the existing manuscript and PNG; no new candidate search was performed.

Changes
- Corrected the script's panel C identifier from DN3635 to DN27036.
- Regenerated both the vector PDF and the 300-dpi PNG from the same script/input.
- Created Figure_9_source_data.tsv with all three panels' displayed data, plus
  panel C adjusted P values, expression status, SignalP/TM and MEROPS fields.
- Changed display spelling from defence to defense for American English.
- Replaced the unavailable external input path with a path beside the script.
- Embedded TrueType fonts in the vector PDF.

Preserved
The input TSV is byte-for-byte identical to the supplied frozen table. All 806
unique candidates, their memberships, statuses and expression statistics are
preserved. Panel A counts remain 86, 28, 518, 14, 3 and 157. Panel B GP/UP counts
remain PR 39/8, jasmonate 5/7, salicylate 2/0, peptidase 40/13 and inhibitor 8/0.
Group memberships overlap; panel B category counts must not be summed.
No DESeq2, annotation, enrichment or sequence-prediction analysis was rerun.
The supplied manuscript was not edited. Its Figure 9 candidate text and caption
require no change to resolve this particular identifier mismatch. Broader
manuscript and figure-layout review remains a separate revision step.

Use
Keep the script and frozen input TSV together. Requires Python 3, matplotlib,
numpy and pandas. Run from any working directory:
  python /path/to/Figure_9_plotting_script.py
Optional:
  python /path/to/Figure_9_plotting_script.py --outdir /path/to/outputs
The script produces the PDF, PNG and source-data TSV in the specified folder,
or beside the script by default. It summarizes existing labels and plots
existing statistics; it does not estimate or change them.

Verification
All ten plotted identifiers match the manuscript selection. Panel C source
statistics preserve the input values. Every plotted candidate is GP_up,
SignalP SP-positive and has has_TMR=FALSE. The PDF contains vector objects and
searchable text with no raster images. Both exported formats were inspected.

Pending next step
RQN 6.3-8.7 is stated in the manuscript, but Supplementary Table S1 and the
sample-level QC report were not supplied in this package. The actual minimum
and affected sample(s) remain unverified; no retention rationale was invented.

Environment used for this plotting export
Python 3.12.14; matplotlib 3.10.8; pandas 2.2.3; numpy 2.3.5
