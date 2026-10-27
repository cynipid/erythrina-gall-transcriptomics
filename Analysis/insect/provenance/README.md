# Retained mapping and assignment evidence

These files document existing runs; no alignment or count assignment was repeated. STAR original records preserve version 2.7.10b, sjdbOverhang 149 and the executed configuration. Original L1/L3/L4 and P3/P4/P5 correspond to manuscript L1/L2/L3 and P1/P2/P3, respectively; see STAR/parsed_log_summary.json.

The featureCounts original-record archive includes the GFF, commands/logs, raw counts, final matrix and metadata. The executed command used `-F GFF -t exon -g locus_tag -s 2`; the retained warning says the GTF parser was used. The audit checks annotation groups/coordinates/union lengths and all 89,796 count cells across 14,966 groups. Assignment percentages use featureCounts-reported Total alignments, which differ from STAR input read pairs. The 152 MB uncompressed annotation is intentionally stored only inside the compressed evidence archive.

To repeat the audit, extract `featureCounts/original_records/EGW_featureCounts_records_comment8.tar.gz` into a new scratch directory and run:

```bash
python Analysis/insect/provenance/featureCounts/audit_featureCounts_records.py --records /path/to/extracted/records --output /path/to/new/audit
```

SOURCE_SHA256SUMS.txt files preserve source identities; original paths inside records are historical provenance rather than portable runtime settings. Audits do not rerun BAM counting or statistical models.
