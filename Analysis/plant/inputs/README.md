# Frozen plant inputs and checkpoint outputs

Fourteen original files supplied on 16 October 2026 are compressed individually with gzip. Decompression restores their exact bytes and original paths relative to the author's `projects` directory. This includes binary RDS files, whose added `.gz` layer must be removed before use. `manifest.tsv` records SHA-256 values before and after compression.

The bundle includes the tximport object, sample metadata, expanded and refined transcriptome-wide mask tables, 292 coding-denominator candidates, earlier expanded-mask results, final result and annotated-result tables, the final fitted DESeq2 object, the annotation donor, GO and KEGG maps, KEGG lookup and frozen expanded-mask defense inventory.

`verification.json` records matches against historical MD5 entries, the 24,734-row final result and status counts, preservation of all result fields after annotation, exclusion of masked genes, and exact agreement of the 16,404 excluded gene IDs with the prepared GEO mask-membership file. These checks read existing files; no statistical models or enrichments were rerun. RDS contents were not loaded in R during this review.

Restore into a new directory using the adjacent `restore_inputs.py`. It decompresses files and verifies checksums; it does not run analysis scripts. The destination must not already exist. For example, from the repository root:

```bash
python Analysis/plant/inputs/restore_inputs.py --destination /path/to/new/EGW_projects
```

This restores the two original project directories beneath that destination. The archived analysis scripts still require path configuration. In particular, several original scripts use `~/projects` directly; `QE_PROJECTS_ROOT` alone does not relocate every dependency. Do not copy over existing analysis projects or execute the historical numbered scripts as an automatic pipeline.

This bundle supplies inputs and checkpoints for the final fit/downstream branch. It does not include every input required by earlier sequence-screening or ranking audits, and it does not by itself recreate the entire provenance history. Audit 39 additionally expects filtering/transition tables and other execution outputs that are not among these 14 files.
