# Erythrina gall transcriptomics

Repository: https://github.com/cynipid/erythrina-gall-transcriptomics

**Preparation draft: figure reproduction package and retained plant analysis scripts. Fourteen frozen plant inputs/checkpoints, plant scripts and execution logs are included. Insect source scripts are also included. Path portability, earlier audit dependencies and additional insect annotation/MEROPS inputs and primary execution records remain incomplete. Selected insect inputs and a retained QC session are now included.**

Associated study: *Parallel transcriptomics reveals insect developmental shifts and gall-associated host transcriptional differences in Quadrastichus erythrinae and Erythrina variegata*.

Author: Serdar Dinç. ORCID: https://orcid.org/0000-0002-8920-2738

## Study design

Samples were collected from one *Erythrina variegata* tree at Academia Sinica, Taipei, Taiwan, on one sampling occasion in November 2020. Insect comparisons comprise three larval and three pupal pools, with ten individuals per pool. Plant comparisons comprise three galled petioles (GP) and three ungalled petioles (UP), with each whole petiole serving as a replicate. The two comparisons were analyzed separately; parallel functional themes do not establish direct molecular interactions.

## Contents

- `Main_Figures/Figure_1` through `Figure_10`: plotting scripts, direct source tables, metadata and required inputs.
- `Analysis/insect/`: original Slurm, DESeq2/enrichment/candidate and prediction-classification scripts, plus blind=TRUE QC sensitivity code.
- `Analysis/plant/`: 40 retained source scripts, a stage map, required-input documentation and checksums. Scripts 38–40 define the final refined-mask branch; scripts 27/28 are dependencies of script 40.
- `rebuild_all.py`: rebuilds the ten main figures using the adjacent files.
- `requirements.txt`: Python plotting dependencies, with versions retained from the figure package where recorded. Pillow is additionally required for the approved Figure 1; its historical version is not specified.
- `docs/source_provenance.tsv`: original local source and SHA-256 for copied files. These paths document provenance, not runtime dependencies.
- `docs/RELEASE_CHECKLIST.md`: remaining tasks before a manuscript release.

The 15 October Figure 1 supersedes the earlier design. Figure 9 retains `TRINITY_DN27036_c0_g1` and includes its source-data TSV. Plotting uses existing estimates; it does not rerun differential expression, enrichment, secretion prediction or MEROPS analyses.

## Reproduce the figures

Use Python 3 in a separate environment. From the repository folder:

```bash
python -m pip install -r requirements.txt
python rebuild_all.py
```

PDF and PNG outputs are written within each figure folder. PCA panels use frozen coordinates; rebuilding them does not refit DESeq2. Figure 9 may regenerate its source-data table from the included frozen candidate overlay.

## Frozen result conventions

Differential-expression criterion: adjusted P < 0.01 and absolute unshrunk log2 fold change >= 1. Insect contrast is pupa versus larva; positive coefficients indicate pupa-up expression. Plant contrast is GP versus UP; positive coefficients indicate GP-up expression.

Final insect results: 4,406 DEGs = 1,909 larva-up + 2,497 pupa-up. Final plant results: 3,530 DEGs = 2,147 GP-up + 1,383 UP-up, after exclusion of the frozen 16,404-gene mask. The plant result universe contains 24,734 genes, including 172 with unavailable adjusted P values. Older plant results must not be substituted.

## Data and code availability status

GEO/SRA deposition is in preparation. No accession has been assigned in this draft. Raw FASTQ files, large sequence files and processed count matrices are not included here. Repository release, commit and public-access details must be finalized before citing code availability in the manuscript.

This package supports figure reproduction from supplied source data. It is not a complete raw-read-to-results pipeline. The final refined-mask plant scripts are included under `Analysis/plant/`; 14 frozen input/checkpoint files and retained final execution logs/session records are included, with checksum verification. Earlier audit dependencies and path portability remain documented limitations. Insect source scripts are included; selected count/QC/prediction inputs and a QC session record are included; additional annotation/MEROPS resources and primary execution records still require consolidation or links to deposited files. Earlier plant scripts are explicitly identified by stage; they do not supersede the final workflow.

## Attribution and licensing

The Figure 1 photograph is reproduced from Lin, Tung and Yang (2021), *Forests* 12, 948, Fig. 1B, under CC BY 4.0; see `Main_Figures/Figure_1/PHOTO_ATTRIBUTION.txt`. https://doi.org/10.3390/f12070948

A software license for the author's code has not yet been selected. The photograph retains its own license. Third-party software retains its respective licenses.
