# Erythrina gall transcriptomics

Code and supporting analysis resources for *Parallel transcriptomics reveals insect developmental shifts and gall-associated host transcriptional differences in Quadrastichus erythrinae and Erythrina variegata*.

Author: Serdar Dinç. [ORCID](https://orcid.org/0000-0002-8920-2738).

Repository: https://github.com/cynipid/erythrina-gall-transcriptomics

## Data deposition

The study and its 12 samples have been approved by NCBI GEO under [GSE351178](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE351178). The deposit includes the raw RNA-seq reads, insect featureCounts table, six plant RSEM gene-level abundance tables, full plant Trinity assembly, the 16,404-gene mask-membership table and sample crosswalk documentation. Raw-read transfer to SRA is pending as of 17 October 2026; SRA links will be available through the GEO record after processing.

GEO records are currently private, with a scheduled public release of 16 October 2027. The release date must be advanced if a publication or preprint citing the accession becomes public earlier. Reviewer access is provided separately to the journal and is not stored in this repository.

## Study design

Samples originated from one *Erythrina variegata* tree at Academia Sinica, Taipei, Taiwan, during one sampling occasion in November 2020. Insect comparisons comprise three larval and three pupal pools, with ten individuals per pool, each from a different galled petiole. Plant comparisons comprise three galled petioles (GP) and three ungalled petioles (UP), with each whole petiole serving as a replicate. The two comparisons were analyzed separately. Parallel functional themes do not establish direct molecular interactions.

## Repository contents

| Location | Contents |
| --- | --- |
| `Main_Figures/Figure_1` through `Figure_10` | Current plotting scripts, direct figure-source data and required inputs. |
| `Analysis/plant/` | Forty original scripts, 14 compressed inputs/checkpoints, final-branch execution logs and version records. |
| `Analysis/insect/` | Twenty-seven original scripts, selected count/QC/prediction inputs and a retained QC session record. |
| `Analysis/annotations_merops/` | Eighty-two annotation/MEROPS source files and a final-status overlay for historical plant priority candidates. |
| `docs/source_provenance.tsv` | Original source paths and identities; these are provenance records, not runtime dependencies. |
| `docs/VERIFICATION.txt` | Recorded packaging and figure-rebuild checks with their limits. |
| `docs/RELEASE_CHECKLIST.md` | Tasks remaining before the manuscript release. |
| `SHA256SUMS.txt` | SHA-256 identities of repository files, excluding this checksum file itself. |

Raw FASTQ reads and the complete plant assembly are deposited through GEO rather than included here. Selected analysis count tables, intermediate objects and annotation resources are included in compressed form; consult each analysis README and manifest for exact contents.

## Reproduce the main figures

Use Python 3 in a separate environment. From the repository root:

```bash
python -m pip install -r requirements.txt
python rebuild_all.py
```

PDF and PNG outputs are written within each figure folder. The plotting requirements are separate from the R, scheduler and external-tool requirements of the archived analyses. PCA panels use frozen coordinates. Figure 9 uses `TRINITY_DN27036_c0_g1` and can regenerate its source table from the supplied candidate overlay. Figure scripts plot existing estimates; they do not rerun differential expression, enrichment, secretion prediction or MEROPS searches.

## Restore archived inputs

From the repository root, choose a new destination for each restoration:

```bash
python Analysis/plant/inputs/restore_inputs.py --destination /path/to/new/EGW_plant_projects
python Analysis/insect/inputs/restore_inputs.py --destination /path/to/new/EGW_insect_projects
python Analysis/annotations_merops/inputs/restore_inputs.py --destination /path/to/new/EGW_annotation_projects
```

The helpers restore original project-relative paths and verify file identities. They do not execute analyses and refuse existing destinations. The added gzip layer must be removed before binary RDS objects are read in R; these helpers perform that restoration.

## Analysis scope and execution limits

The package supports reproduction of the main figures from supplied inputs and preserves selected analysis sources, inputs, results and execution evidence. It is not a turnkey raw-read-to-results workflow. Archived scripts retain original Windows/WSL/TRUBA path assumptions and require configuration, external software and some inputs not included here. A retained container-build job depends on a local sandbox that is not supplied. Third-party model weights and the MEROPS database are not distributed.

Plant scripts 38–40 document the final refined-mask branch; script 40 adapts scripts 27/28 in memory. Running 27/28 directly targets an earlier branch. Do not run the full numbered script collection as an automatic pipeline. Earlier analyses and historical MEROPS tables do not supersede final manuscript results. The MEROPS README explains the final-status overlay: 19 of 20 historical plant priority candidates remain GP-up; `TRINITY_DN12519_c0_g1` is Non-DE in the final fit.

Original source scripts and data were retained unchanged during documentation updates. Packaging checks did not rerun statistical models, alignments, enrichment tests or prediction/search tools. The retained R-session files document specific historical sessions rather than a complete environment lockfile.

## Final result conventions

Differential expression requires adjusted P < 0.01 and absolute unshrunk log2 fold change >= 1. Insect coefficients compare pupa versus larva; positive values indicate pupa-up expression. Plant coefficients compare GP versus UP; positive values indicate GP-up expression.

- Insect: 4,406 DEGs, comprising 1,909 larva-up and 2,497 pupa-up genes; 13,919 retained genes, including four unavailable adjusted P values.
- Plant: 3,530 DEGs, comprising 2,147 GP-up and 1,383 UP-up genes after exclusion of the frozen 16,404-gene mask; 24,734 retained genes, including 172 unavailable adjusted P values.

## Attribution and licensing

The Figure 1 photograph is reproduced from Lin, Tung and Yang (2021), *Forests* **12**, 948, Fig. 1B, under CC BY 4.0. See `Main_Figures/Figure_1/PHOTO_ATTRIBUTION.txt` and https://doi.org/10.3390/f12070948.

The author's original analysis scripts, plotting scripts and utility code are licensed under the MIT License; see `LICENSE`. Copyright (c) 2026 Serdar Dinç. The MIT License applies to original code and its associated documentation only. It does not relicense research data, third-party software, database-derived annotations, reference resources or the Figure 1 photograph. Their applicable terms and attribution remain in force.
