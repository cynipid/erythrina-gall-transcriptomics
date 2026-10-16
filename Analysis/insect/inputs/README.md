# Insect inputs and retained QC session

Five files from the original `qe_genome_guided` project are compressed individually. The manifest records their original project-relative paths and SHA-256 values. The additional gzip layer must be removed before R reads the VST RDS file. To restore exact bytes into a new directory, from the repository root:

```bash
python Analysis/insect/inputs/restore_inputs.py --destination /path/to/new/EGW_insect_projects
```

The helper does not run analysis scripts and refuses an existing destination. Original archived analysis scripts still need path configuration; see the parent README.

Read-only checks found 14,966 count rows and six columns matching the sample metadata. All 89,796 integer count values exactly match the prepared GEO featureCounts source after explicit sample-column mapping. The total-count >=10 filter retains 13,919 gene IDs, exactly matching the Figure 3 source universe. The 13,455-row prediction merge table has the SHA-256 expected by the custom classification script. No DESeq2 fit, transformation, secretion prediction or classification was rerun.

The retained `R_session_info.txt` documents the reviewer-requested blind=TRUE QC environment: R 4.6.1, DESeq2 1.52.0, ggplot2 4.0.3 and SummarizedExperiment 1.42.0, on Windows 11 x64. It does not establish historical STAR, SAMtools, featureCounts or prediction-tool versions. The VST RDS is retained byte-for-byte but was not loaded into R during this review.

This is a selected input bundle. Annotation mappings, candidate source data, retained primary result/environment records and MEROPS resources are not supplied by these five files. The repository is not yet a complete insect analysis execution bundle.
