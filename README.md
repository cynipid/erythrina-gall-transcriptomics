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
| `Main_Figures/Figure_1` through `Figure_6` | Final PDF/PNG artwork, captions, portable plotting scripts and frozen source inputs for all six main figures. Figure 1 retains its full-resolution photo input. |
| `Supplementary_Figures/Figure_S1` through `Figure_S8` | Final supplementary PDF/PNG files, captions, plotting scripts and frozen source inputs. |
| `Main_Figures/Main_Figure_Captions.docx` | Editable main captions preserving scientific italics and subscripts. |
| `Main_Figures/` and `Supplementary_Figures/` caption/preview builders | `captions.json` is the canonical formatted caption source; `build_captions.py` writes individual captions, combined TXT and editable DOCX, and `build_preview.py` builds each combined PDF. |
| `docs/Figure_Numbering_Map.tsv` | Previous-to-final numbering map for all fourteen figures. |
| `Analysis/plant/` | Forty original scripts, 14 compressed inputs/checkpoints, final-branch execution logs and version records. |
| `Analysis/insect/` | Twenty-seven original scripts, selected count/QC/prediction inputs and a retained QC session record. |
| `Analysis/annotations_merops/` | Eighty-two annotation/MEROPS source files and a final-status overlay for historical plant priority candidates. |
| `Analysis/insect/provenance/` | Retained STAR and featureCounts records and the corresponding mapping/count-assignment audits. |
| `Analysis/plant/sample_expression/` | A 32-cluster sample-level expression extract, portable extraction script, frozen normalized-count/statistics inputs and validation records. |
| `Analysis/supplementary_tables/` | Revised Supplementary Tables S2 and S13, with their scope and source links documented. |
| `docs/source_provenance.tsv` | Historical original-snapshot source paths and identities; some recorded repository paths use the former ten-main-figure organization. See `docs/README.md`. |
| `docs/VERIFICATION.txt` | Recorded packaging and figure-rebuild checks with their limits. |
| `docs/RELEASE_CHECKLIST.md` | Completed synchronization and remaining submission/archive follow-up. |
| `SHA256SUMS.txt` | SHA-256 identities of repository files, excluding this checksum file itself. |

Raw FASTQ reads and the complete plant assembly are deposited through GEO rather than included here. Selected analysis count tables, intermediate objects and annotation resources are included in compressed form; consult each analysis README and manifest for exact contents.

## Final figure organization

The final manuscript has six main figures and eight supplementary figures. The insect and plant PCA and volcano plots are supplementary. Every numbered folder contains its current caption, a portable `plot.py` and frozen source data. Final PDF and PNG artwork is included for all six main figures and all eight supplementary figures, together with both combined preview PDFs. The same publication masters are supplied in the final figure ZIP packages. The four special-upload master-file paths are listed in `docs/MASTER_FILE_UPLOAD_LIST.tsv`; this crosswalk is not the complete repository inventory. `Figure_Index.tsv` in each figure collection lists the exact files. The earlier ten-main-figure layout is mapped in `docs/Figure_Numbering_Map.tsv`.

The 15 October Figure 1 revision is retained. The 26 October update enlarges the complete identifiers in main Figures 5 and 6 while preserving candidate selection, values, classes, colors, order and axis ranges. Figure 5 has three vertically stacked panels at 180 × 203 mm with 9 pt identifiers; Figure 6 retains four panels at 180 × 240 mm with 10 pt identifiers in panel c. Both combined previews have nine pages. Main Figures 5 and 6 each have a separate following caption page so their artwork remains large; their identifier fonts are approximately 8.81 and 9.79 pt in the main preview. Supplementary correlation plots S2 and S6 retain their portable layouts with unchanged coefficients and sample labels. S2 includes its recovered full-precision correlation matrix, frozen transformed input, sample crosswalk and documented optional extraction script. Statistical models were not refitted.

The existing v1.0.0 and v1.1.0 releases remain earlier frozen snapshots. This 26 October synchronization updates the current figure display, approved captions and supporting provenance without refitting statistical models. Cite the exact resulting commit using its full commit identifier (the GitHub `/tree/<commit-SHA>` URL), or an explicitly created later archived release. Do not cite an earlier release as containing these revisions.

## Reproduce the figures

Use Python 3 in a separate environment. From the repository root:

```bash
python -m pip install -r requirements.txt
python rebuild_all.py --outdir /path/to/new/EGW_rebuilt_figures
```

The command rebuilds all fourteen individual figures. Use `--group main` or `--group supplementary` for one collection. Without `--outdir`, outputs are written beside the supplied source files. Each collection also has its own `rebuild_all.py` and `build_preview.py`. The preview builder uses the figure PDFs and canonical formatted captions in that collection. To regenerate editable and plain-text captions after an approved change to `captions.json`, run `python Main_Figures/build_captions.py` or `python Supplementary_Figures/build_captions.py`. Caption generation requires `python-docx==1.2.0`, included in the requirements.

The 17 October checks executed all fourteen plot scripts in isolated output directories against that snapshot's production artwork. Current export verification is recorded in each collection's `Export_Audit.json`: all seven main source TSVs remain byte-identical to the preceding approved package, main Figures 1–4 remain unchanged, and all 55 non-caption supplementary figure-folder files remain byte-identical. Figure 5/6 artwork reflects the authorized readability changes described above. The historical `Reproduction_Audit.json` and earlier records are retained as dated evidence, rather than as pixel-identity claims about newly revised artwork. PDF metadata may vary between runs. PCA panels use frozen coordinates; correlation panels plot supplied full-precision matrices. Main Figure 5 contains `TRINITY_DN27036_c0_g1` and regenerates its source summary from the supplied candidate overlay. Ordinary plotting does not require RDS conversion. Optional S2 source extraction requires `rdata==1.1.0`.

The plotting requirements are separate from the R, scheduler and external-tool requirements of the archived analyses. Figure scripts use existing estimates and classifications; they do not rerun differential-expression models, enrichment, secretion prediction or MEROPS searches.

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

The historical analysis archive is retained. New STAR/featureCounts evidence and the plant sample-level expression extract are documented in their dedicated subfolders; the revised S2 and S13 workbooks are included separately. S13 adds a sample-level worksheet for 32 unique plant clusters from frozen final normalized counts and statistics, without altering the original 806-candidate inventory. These additions do not constitute a new analysis fit. The 26 October figure update changes display scripts and caption/preview builders, not frozen plotted source values. Packaging and extraction checks did not rerun statistical models, alignments, enrichment tests or prediction/search tools. The retained R-session files document specific historical sessions rather than a complete environment lockfile.

## Final result conventions

Differential-expression classification requires a BH-adjusted *P* value < 0.01 together with an absolute estimated, unshrunk log2 fold change >= 1. The fold-change criterion is a filter on the fitted point estimate, not a formal test against a nonzero effect-size threshold. Insect coefficients compare pupa versus larva; positive values indicate pupa-up expression. Plant coefficients compare GP versus UP; positive values indicate GP-up expression.

- Insect: 4,406 DEGs, comprising 1,909 larva-up and 2,497 pupa-up genes; 13,919 retained genes, including four unavailable adjusted P values.
- Plant: 3,530 DEGs, comprising 2,147 GP-up and 1,383 UP-up genes after exclusion of the frozen 16,404-gene mask; 24,734 retained genes, including 172 unavailable adjusted P values.

## Attribution and licensing

The Figure 1 photograph is reproduced from Lin, Tung and Yang (2021), *Forests* **12**, 948, Fig. 1B, under CC BY 4.0. See `Main_Figures/Figure_1/PHOTO_ATTRIBUTION.txt` and https://doi.org/10.3390/f12070948.

The author's original analysis scripts, plotting scripts and utility code are licensed under the MIT License; see `LICENSE`. Copyright (c) 2026 Serdar Dinç. The MIT License applies to original code and its associated documentation only. It does not relicense research data, third-party software, database-derived annotations, reference resources or the Figure 1 photograph. Their applicable terms and attribution remain in force.

