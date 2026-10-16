# Annotation resources and MEROPS source snapshot

This addition preserves 82 source files selected from the author's existing projects. Eight MEROPS scripts remain readable as source under `inputs/projects/qe_merops_screen/90_manuscript_release/04_current_scripts/`; all other source files are individually gzip-compressed. `inputs/manifest.tsv` records original paths, exact sizes and SHA-256 values. No original source contents were changed.

## Important version distinction

The MEROPS source is the frozen 12 September 2026 curated-v2 snapshot. Its README and plant result tables predate the final 16,404-gene mask. Their historical headline counts and plant expression labels must not replace final manuscript values.

The supplied 20 historical plant manuscript-priority candidates were joined to the existing final refined-mask result table, without refitting or searching. Nineteen remain GP-up and exactly match Figure 10C. `TRINITY_DN12519_c0_g1` is Non-DE in the final result. The derived `plant_priority_final_status_overlay.tsv` records all 20 with final coefficients, standard errors and adjusted P values. Other historical candidate-status rows are retained unchanged and must not be interpreted as final-status overlays.

The historical insect inhibitor-family homologue set contains 149 unique genes. Its nine-gene larval shortlist exactly matches Figure 10B. MEROPS homology and secretion predictions do not establish biochemical inhibition, host exposure, plant targets or direct molecular interactions. The plant search was restricted to its original GP-up classical-secretory-pathway query set; this archive does not expand that search or reselect query proteins.

## Annotation resources

The insect files include original eggNOG output, parsed mappings and earlier full/named enrichment tables. The plant files include original TransDecoder-protein eggNOG output, the representative-per-gene subset and processed gene annotations. Headers of all three supplied `.annotations` files record eggNOG-mapper 2.1.13; this header evidence does not itself establish an eggNOG database version. Original command lines and timestamps remain in the decompressed files.

## Restore exact original files

From the repository root:

```bash
python Analysis/annotations_merops/inputs/restore_inputs.py --destination /path/to/new/EGW_annotation_projects
```

Choose a new destination. The helper restores and checks exact original bytes, including readable source scripts. It runs no analyses. The restored tree uses original project-relative paths; archived scripts retain original path assumptions and are not yet portable. The MEROPS database itself is not included; its retained metadata/checksums are supplied under the restored `08_database_metadata` directory. Nineteen files referenced by the original full-release manifest were not selected for this addition; those absences are recorded in `inputs/verification.json`, not treated as checksum failures.

## Verification scope

All 57 selected files referenced by the MEROPS release manifest matched its SHA-256 values. All 82 stored files were checked against their originals after compression/copying. Four Python scripts passed syntax parsing; R syntax was not checked because R was unavailable. This review compared candidate identities and joined existing final results only. It did not rerun MEROPS, DIAMOND, eggNOG, enrichment, secretion prediction or DESeq2. Historical scripts require additional input tables and path configuration; this is not a complete executable pipeline.
