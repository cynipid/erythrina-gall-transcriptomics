# Plant analysis scripts: retained source archive

All 40 source files are preserved byte-for-byte from `EGW_analysis_scripts_for_review.tar.gz`, supplied by the author on 16 October 2026. The inventory records roles and checksums. This is a code archive with documented dependencies, not a standalone execution bundle.

## Final manuscript branch

| Script | Role |
| --- | --- |
| 37 | Reviews nucleotide coverage using coding-aware denominators; identifies the 292 additions to the expanded mask. |
| 38 | Unions the 16,112-gene expanded mask and 292 additional genes; removes all 16,404 genes from tximport matrices before DESeqDataSetFromTximport, filtering and fitting. |
| 39 | Independently checks final result integrity, mask exclusion and transitions against the expanded-mask result. |
| 40 | Joins frozen annotations, refreshes QC and enrichment, curates GO labels and overlays the frozen 806-gene defense inventory. |
| 27, 28 | Required by script 40: their source is adapted in memory to the final refined-mask inputs and expected counts. Running these two directly targets the earlier expanded-mask branch. |

Do not execute scripts 01–40 as a new automatic pipeline. Earlier scripts document successive review stages and optional discovery rankings; their values do not replace the final manuscript result. No scripts were executed for this repository packaging task.

Final result controls: 118,082 original genes; 16,404 masked; 101,678 genes before the low-count filter; 24,734 after requiring counts >=10 in >=3 samples; 2,147 GP-up, 1,383 UP-up, 21,032 Non-DE and 172 unavailable adjusted P values. Model `~ condition`; contrast GP versus UP; padj <0.01 and |unshrunk log2FC| >=1. Retained source metadata use `SG` = GP and `S` = UP. Historical variable names and data columns were preserved.

## Required input data

These scripts refer to the author's original adjacent projects `qe_plant_denovo` and `qe_cross_species_discovery`. Those project inputs are not included in this script-only addition. In particular, script 38 requires:

- `qe_plant_denovo/07_deseq2/plant_gene_txi.rds` (counts, abundance and length matrices).
- `qe_plant_denovo/01_input/expression/plant_deseq2_samples.tsv`.
- `qe_cross_species_discovery/04_cross_species_integration/results/plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv`.
- `qe_cross_species_discovery/04_cross_species_integration/results/plant_coding_denominator_additional_sequence_review_candidates.tsv`.
- `qe_cross_species_discovery/05_host_focused_DESeq2_sensitivity/results/plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_all.tsv`.

Script 40 additionally requires the final result and fitted DESeq2 object produced by script 38, the transcriptome-wide refined-mask table, the original eggNOG annotation donor, the frozen expanded-mask defense inventory, and GO/KEGG gene mappings and pathway-name lookup used by script 28. All must be reconciled with deposited data or supplied as supporting inputs before claiming a complete reproducible analysis package. The 16,404-row GEO mask-membership file is not a drop-in replacement for the 118,082-row transcriptome-wide mask table expected by these scripts.

## Paths, environments and verification limits

Scripts 38 and 40 support `QE_PROJECTS_ROOT`; several dependencies and Python audits use `~/projects` or the original WSL path directly. Setting that environment variable alone does not relocate every script, including script 28. Portability changes remain pending; no source paths or scientific logic were silently rewritten.

Dependencies visible in the final branch include R, DESeq2, tximport, SummarizedExperiment, ggplot2, AnnotationDbi and GO.db; Python audits use Python 3. This list is not a locked environment for every historical script. Historical `sessionInfo()` records and retained audit logs should be supplied before completing environment documentation.

Static inspection confirmed the final masking order, model, count filter, contrast and classification rules. All 31 Python files passed syntax parsing. R was not available for syntax parsing here; no analysis was rerun and no new output counts were calculated. Assertions in source code describe expected controls, not evidence of execution in this packaging review.
