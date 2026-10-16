# Insect analysis source archive

Twenty-seven source files are preserved byte-for-byte from `EGW_insect_scripts_for_review.tar.gz`, supplied by the author on 16 October 2026. No scripts were executed or jobs submitted during packaging. `script_inventory.tsv` records identities and syntax-check limits.

## Contents and roles

- `90_manuscript_release/06_current_scripts/truba_slurm`: five retained Slurm job scripts for input checks, container construction, STAR indexing/alignment and featureCounts strand comparisons.
- `90_manuscript_release/06_current_scripts/deseq2_enrichment_candidates`: twenty R scripts covering differential expression, annotation, GO/KEGG enrichment, candidate inventories and historical plotting.
- `90_manuscript_release/06_current_scripts/signalp_deeptmhmm`: one custom script classifying existing SignalP/DeepTMHMM prediction outputs.
- `07_deseq2/insect_blind_TRUE_QC_sensitivity_RStudio.R`: additional reviewer-requested QC sensitivity script.

No installed SignalP package, model weights or third-party source distribution is included. The custom classification script operates on frozen prediction tables; it does not run SignalP or DeepTMHMM predictions.

## Final analysis conventions

The primary script `02_run_deseq2_primary_FDR001_LFC1.R` constructs a DESeq2 dataset with `~ stage`, retains genes with total counts >=10 across six libraries, fits DESeq2 and extracts pupa versus larva with `alpha=0.01`. Reporting requires padj <0.01 and |unshrunk log2FoldChange| >=1. Positive values indicate pupa-up expression. The final manuscript result is 4,406 DEGs: 1,909 larva-up and 2,497 pupa-up; 13,919 genes are retained, including four unavailable adjusted P values. These output controls have not been recalculated in this script-only review.

The archived STAR index job specifies `sjdbOverhang=149` and `genomeSAindexNbases=13`. The alignment job supplies `readFilesCommand=zcat` and `outSAMtype=BAM Unsorted`; it does not explicitly specify STAR two-pass mode. The featureCounts job compares strand settings 0, 1 and 2, using `-p --countReadPairs -B -C -F GFF -t exon -g locus_tag`; the manuscript uses reverse-stranded counts (`-s 2`). Script presence alone does not establish installed tool versions or historical execution.

The custom secretory-pathway classifier specifies 13,455 proteins and expects 2,039 classical-secretory-pathway candidates; downstream-TMR and discordant cases are excluded from its primary group. This prediction classification does not establish secretion or effector activity experimentally.

## Sample identifiers

| Original sample | Manuscript label |
| --- | --- |
| L1 | L1 |
| L3 | L2 |
| L4 | L3 |
| P3 | P1 |
| P4 | P2 |
| P5 | P3 |

Original filenames and metadata identifiers remain unchanged. The sensitivity script records this display mapping and compares blind=TRUE QC with the retained primary blind=FALSE transformation; it does not refit the primary differential-expression model.

## Inputs, paths and pending consolidation

The primary DESeq2 script requires `qe_genome_guided/06_counts/gene_counts_matrix.tsv` and `sample_metadata.tsv`. The sensitivity script additionally requires `08_results_primary/insect_gene_vsd.rds`. The prediction classifier requires `11_signalp_tmhmm/results/signalp6_deeptmhmm_genomewide_13455.tsv`. Annotation and enrichment scripts require their named annotation/mapping inputs. These four named inputs and a retained QC R-session record are now provided under `inputs/`; see its README and verification report. Additional annotation/enrichment resources are supplied separately under `Analysis/annotations_merops/`; some runtime dependencies and primary-analysis execution records remain outside the package.

Slurm jobs retain the original `/arf/scratch/sdinc/qe_rnaseq_truba` paths and scheduler settings. Several R scripts use the author's Windows-to-WSL project path. The container build job refers to a local sandbox not supplied here. These are original execution records, not portable jobs ready to submit from a repository checkout. No versions have been inferred from filenames or container names.

Current manuscript figure reproduction is provided separately under `Main_Figures`; those scripts supersede historical plotting layouts in this source archive. MEROPS and annotation resources are documented separately under `Analysis/annotations_merops/`. This addition does not constitute a complete raw-read-to-results workflow.

## Packaging verification

The one Python file passed syntax parsing. All five Slurm files passed `bash -n` syntax checks. R was unavailable for syntax parsing. Source files were copied unchanged and SHA-256 values recorded. No DESeq2, enrichment, alignment, featureCounts, SignalP, DeepTMHMM or MEROPS analysis was run.

## Supplied input checks

All six count columns match the prepared GEO source exactly, metadata agree, and the 13,919-gene low-count-filtered universe matches the main Figure 3 source. The objective prediction table matches the classifier input SHA-256. `software_versions.tsv` explicitly limits the recovered R versions to the retained QC session.
