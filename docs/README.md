# Documentation scope

Current figure and caption export checks are recorded in `Main_Figures/Export_Audit.json`, `Supplementary_Figures/Export_Audit.json` and `VERIFICATION.txt`. Earlier reproduction and production audits are preserved with their original dates; their artwork comparisons and file counts apply to those snapshots.

`source_provenance.tsv` preserves original source paths and identities from the earlier ten-main-figure package. Its `repository_file` entries are historical identifiers, and some no longer name current paths after figure renumbering. They are not runtime dependencies or a current file inventory. Use `Figure_Numbering_Map.tsv`, the collection `Figure_Index.tsv` files and the current root `SHA256SUMS.txt` to locate and verify the synchronized files.

`MASTER_FILE_UPLOAD_LIST.tsv` is a four-artifact special-upload crosswalk; it is not a complete file list. Its ZIP filenames retain the original package naming convention, while current package contents and checks are identified by their manifests and export-audit dates.

The 26 October supporting additions are documented at:

- [STAR provenance](../Analysis/insect/provenance/STAR/)
- [featureCounts provenance](../Analysis/insect/provenance/featureCounts/)
- [Plant sample-level expression](../Analysis/plant/sample_expression/)
- [Revised Supplementary Tables S2 and S13](../Analysis/supplementary_tables/)

These records document retained analyses and a frozen-data expression extraction; they do not imply new alignment, count-assignment or statistical-model runs. Confidential GEO reviewer-access information is distributed separately and is not part of this public archive.

