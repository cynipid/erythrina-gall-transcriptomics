#!/usr/bin/env Rscript

# Final refined-mask host-focused GP-versus-UP DESeq2 analysis.
#
# This script freezes the 16,112-gene expanded nucleotide mask and adds the
# 292 coding-denominator-supported genes identified by the completed,
# transcriptome-wide audit. It performs one final DESeq2 refit and writes only
# new files containing "refined_mask". No earlier result is overwritten.

options(stringsAsFactors = FALSE)

required_packages <- c("DESeq2", "tximport", "SummarizedExperiment")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(SummarizedExperiment)
})

projects_root <- Sys.getenv("QE_PROJECTS_ROOT", unset = "")
if (!nzchar(projects_root)) {
  projects_root <- if (.Platform$OS.type == "windows") {
    "//wsl.localhost/Ubuntu/home/dodi/projects"
  } else {
    file.path(Sys.getenv("HOME"), "projects")
  }
}

discovery_root <- file.path(projects_root, "qe_cross_species_discovery")
plant_root <- file.path(projects_root, "qe_plant_denovo")
integration_dir <- file.path(discovery_root, "04_cross_species_integration", "results")
branch_root <- file.path(discovery_root, "05_host_focused_DESeq2_sensitivity")
results_dir <- file.path(branch_root, "results")
logs_dir <- file.path(branch_root, "logs")
control_dir <- file.path(discovery_root, "00_project_control")

input_txi <- file.path(plant_root, "07_deseq2", "plant_gene_txi.rds")
input_metadata <- file.path(
  plant_root, "01_input", "expression", "plant_deseq2_samples.tsv"
)
input_expanded_mask <- file.path(
  integration_dir,
  "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv"
)
input_coding_candidates <- file.path(
  integration_dir,
  "plant_coding_denominator_additional_sequence_review_candidates.tsv"
)
input_expanded_results <- file.path(
  results_dir,
  "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_all.tsv"
)

outputs <- list(
  refined_mask = file.path(
    integration_dir,
    "plant_refined_coding_denominator_sequence_review_mask.tsv"
  ),
  all_results = file.path(
    results_dir,
    "plant_gene_DESeq2_GP_vs_UP_host_focused_refined_mask_all.tsv"
  ),
  normalized_counts = file.path(
    results_dir,
    "plant_gene_normalized_counts_host_focused_refined_mask.tsv"
  ),
  filtering_summary = file.path(
    results_dir,
    "plant_host_focused_refined_mask_DESeq2_filtering_summary.tsv"
  ),
  mask_exclusion_audit = file.path(
    results_dir,
    "plant_refined_coding_denominator_mask_exclusion_audit.tsv"
  ),
  expanded_comparison = file.path(
    results_dir,
    "plant_host_focused_expanded_vs_refined_mask_DESeq2_comparison.tsv"
  ),
  dds = file.path(
    results_dir, "plant_gene_dds_host_focused_refined_mask.rds"
  ),
  parameters = file.path(
    control_dir,
    "38_refined_mask_host_focused_plant_DESeq2_parameters.tsv"
  ),
  input_identity = file.path(
    logs_dir, "38_host_focused_refined_mask_DESeq2_input_identity.tsv"
  ),
  session = file.path(
    logs_dir, "38_host_focused_refined_mask_DESeq2_sessionInfo.txt"
  ),
  audit = file.path(
    logs_dir, "38_host_focused_refined_mask_DESeq2_audit.txt"
  )
)

audit_title <- "REFINED-MASK HOST-FOCUSED PLANT DESEQ2 FINAL ANALYSIS"

error_exit <- function(message) {
  cat(
    "===== ", audit_title, " =====\n",
    "ERROR: ", message, "\n",
    audit_title, ": CHECK REQUIRED\n",
    sep = "", file = stderr()
  )
  quit(status = 1L)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(control_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(outputs$audit)) {
  previous_audit <- readLines(outputs$audit, warn = FALSE)
  if (length(previous_audit) > 0L &&
      tail(previous_audit, 1L) == paste0(audit_title, ": PASS")) {
    error_exit(paste(
      "A completed validated audit already exists:", outputs$audit,
      "No outputs were changed."
    ))
  }
}

required_inputs <- c(
  input_txi, input_metadata, input_expanded_mask,
  input_coding_candidates, input_expanded_results
)
for (path in required_inputs) {
  if (!file.exists(path)) error_exit(paste("Missing required input:", path))
}

read_tsv <- function(path, ...) {
  read.delim(
    path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("NA"), ...
  )
}

write_tsv <- function(object, path, row.names = FALSE) {
  write.table(
    object, file = path, sep = "\t", quote = FALSE,
    row.names = row.names, col.names = TRUE, na = "NA"
  )
}

input_md5_before <- unname(tools::md5sum(required_inputs))
names(input_md5_before) <- basename(required_inputs)

txi <- readRDS(input_txi)
if (!all(c("abundance", "counts", "length") %in% names(txi))) {
  error_exit("plant_gene_txi.rds does not contain abundance, counts and length")
}
if (is.null(rownames(txi$counts)) || is.null(colnames(txi$counts))) {
  error_exit("tximport counts lack gene or sample identifiers")
}
if (anyDuplicated(rownames(txi$counts)) > 0L) {
  error_exit("tximport counts contain duplicate gene IDs")
}
original_gene_ids <- rownames(txi$counts)
if (length(original_gene_ids) != 118082L) {
  error_exit(paste(
    "Expected 118082 genes in the frozen tximport object; found",
    length(original_gene_ids)
  ))
}

metadata <- read_tsv(input_metadata)
if (!all(c("sample_id", "condition") %in% names(metadata))) {
  error_exit("Sample metadata must contain sample_id and condition")
}
if (anyDuplicated(metadata$sample_id) > 0L) {
  error_exit("Duplicate sample_id in metadata")
}
rownames(metadata) <- metadata$sample_id
if (!identical(colnames(txi$counts), rownames(metadata))) {
  error_exit("tximport sample columns do not exactly match metadata sample order")
}
if (!setequal(unique(metadata$condition), c("S", "SG"))) {
  error_exit("Expected frozen legacy metadata conditions S and SG")
}
# Legacy labels: SG = GP and S = UP.
metadata$condition <- factor(metadata$condition, levels = c("S", "SG"))

mask <- read_tsv(input_expanded_mask)
needed_mask_columns <- c("gene_id", "proposed_expanded_mask_status")
if (!all(needed_mask_columns %in% names(mask))) {
  error_exit(paste(
    "Expanded mask lacks required columns:",
    paste(setdiff(needed_mask_columns, names(mask)), collapse = ", ")
  ))
}
if (nrow(mask) != 118082L || anyDuplicated(mask$gene_id) > 0L) {
  error_exit("Expanded mask must contain 118082 unique Trinity gene IDs")
}
if (!setequal(mask$gene_id, original_gene_ids)) {
  error_exit("Expanded mask and tximport gene universes are not identical")
}

expanded_exclusion_statuses <- c(
  "review_exclude_existing_confirmed_multilevel_sequence",
  "review_exclude_additional_strong_nucleotide_similarity"
)
expanded_mask_genes <- unique(mask$gene_id[
  mask$proposed_expanded_mask_status %in% expanded_exclusion_statuses
])
if (length(expanded_mask_genes) != 16112L) {
  error_exit(paste(
    "Expected 16112 genes in the frozen expanded mask; found",
    length(expanded_mask_genes)
  ))
}

coding <- read_tsv(input_coding_candidates)
needed_coding_columns <- c(
  "gene_id", "current_expanded_mask_DESeq2_status",
  "coding_denominator_review_rules", "proposed_review_status"
)
if (!all(needed_coding_columns %in% names(coding))) {
  error_exit(paste(
    "Coding-denominator candidate table lacks required columns:",
    paste(setdiff(needed_coding_columns, names(coding)), collapse = ", ")
  ))
}
if (nrow(coding) != 292L || anyDuplicated(coding$gene_id) > 0L) {
  error_exit(paste(
    "Coding-denominator candidate table must contain 292 unique genes; found",
    nrow(coding)
  ))
}
if (!all(coding$gene_id %in% original_gene_ids)) {
  error_exit("One or more coding-denominator candidates are absent from tximport")
}
if (length(intersect(coding$gene_id, expanded_mask_genes)) > 0L) {
  error_exit("Coding-denominator candidates unexpectedly overlap the expanded mask")
}
if (!all(coding$proposed_review_status ==
         "review_add_coding_denominator_supported_sequence")) {
  error_exit("Unexpected proposed_review_status in coding-denominator candidates")
}
if (any(is.na(coding$coding_denominator_review_rules)) ||
    any(!nzchar(coding$coding_denominator_review_rules))) {
  error_exit("Every coding-denominator candidate must have at least one rule")
}

expected_candidate_statuses <- c(
  GP_up = 192L,
  UP_up = 0L,
  non_DE = 3L,
  NA_padj_not_testable = 2L,
  not_in_expanded_mask_host_focused_tested_universe = 95L
)
observed_candidate_statuses <- table(coding$current_expanded_mask_DESeq2_status)
for (status in names(expected_candidate_statuses)) {
  observed <- if (status %in% names(observed_candidate_statuses)) {
    unname(observed_candidate_statuses[status])
  } else {
    0L
  }
  if (observed != expected_candidate_statuses[status]) {
    error_exit(paste(
      "Coding-candidate status mismatch for", status,
      ": observed", observed, "expected", expected_candidate_statuses[status]
    ))
  }
}
unexpected_candidate_statuses <- setdiff(
  names(observed_candidate_statuses), names(expected_candidate_statuses)
)
if (length(unexpected_candidate_statuses) > 0L) {
  error_exit(paste(
    "Unexpected coding-candidate DESeq2 status:",
    paste(unexpected_candidate_statuses, collapse = ", ")
  ))
}

expanded_results <- read_tsv(input_expanded_results)
needed_expanded_columns <- c(
  "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
  "pvalue", "padj", "DESeq2_status_host_focused_expanded_mask"
)
if (!all(needed_expanded_columns %in% names(expanded_results))) {
  error_exit("Expanded-mask DESeq2 table lacks required columns")
}
if (nrow(expanded_results) != 24931L ||
    anyDuplicated(expanded_results$gene_id) > 0L) {
  error_exit("Expanded-mask DESeq2 table must contain 24931 unique genes")
}
expected_expanded_counts <- c(
  GP_up = 2349L,
  UP_up = 1386L,
  non_DE = 21022L,
  NA_padj_not_testable = 174L
)
observed_expanded_counts <- table(
  expanded_results$DESeq2_status_host_focused_expanded_mask
)
for (status in names(expected_expanded_counts)) {
  observed <- if (status %in% names(observed_expanded_counts)) {
    unname(observed_expanded_counts[status])
  } else {
    0L
  }
  if (observed != expected_expanded_counts[status]) {
    error_exit(paste(
      "Expanded-mask status mismatch for", status,
      ": observed", observed, "expected", expected_expanded_counts[status]
    ))
  }
}
if (length(intersect(expanded_results$gene_id, expanded_mask_genes)) > 0L) {
  error_exit("Expanded-mask results unexpectedly contain an expanded-mask gene")
}

coding_tested <- coding[coding$gene_id %in% expanded_results$gene_id, , drop = FALSE]
if (nrow(coding_tested) != 197L) {
  error_exit(paste(
    "Expected 197 coding-denominator candidates in expanded results; found",
    nrow(coding_tested)
  ))
}
for (gene_id in coding_tested$gene_id) {
  expected_status <- coding$current_expanded_mask_DESeq2_status[
    match(gene_id, coding$gene_id)
  ]
  observed_status <- expanded_results$DESeq2_status_host_focused_expanded_mask[
    match(gene_id, expanded_results$gene_id)
  ]
  if (!identical(expected_status, observed_status)) {
    error_exit(paste(
      "Candidate/expanded-result status disagreement for", gene_id
    ))
  }
}

coding_genes <- coding$gene_id
refined_mask_genes <- union(expanded_mask_genes, coding_genes)
if (length(refined_mask_genes) != 16404L) {
  error_exit(paste(
    "Expected 16404 genes in the refined mask; found",
    length(refined_mask_genes)
  ))
}

# Freeze an explicit transcriptome-wide refined-mask table.
refined_mask <- mask
refined_mask$additional_coding_denominator_mask_overlap <- ifelse(
  refined_mask$gene_id %in% coding_genes, "TRUE", "FALSE"
)
refined_mask$refined_mask_status <- ifelse(
  refined_mask$gene_id %in% expanded_mask_genes,
  "review_exclude_frozen_expanded_nucleotide_mask",
  ifelse(
    refined_mask$gene_id %in% coding_genes,
    "review_exclude_additional_coding_denominator_supported_sequence",
    "retain_outside_refined_sequence_review_mask"
  )
)
rule_lookup <- setNames(coding$coding_denominator_review_rules, coding$gene_id)
refined_mask$coding_denominator_review_rules <- unname(
  rule_lookup[refined_mask$gene_id]
)
refined_mask$coding_denominator_review_rules[
  is.na(refined_mask$coding_denominator_review_rules)
] <- ""
write_tsv(refined_mask, outputs$refined_mask)

# Remove the complete refined mask before DESeq2 construction and fitting.
keep_gene_ids <- original_gene_ids[!original_gene_ids %in% refined_mask_genes]
txi_refined <- txi
for (component_name in names(txi_refined)) {
  component <- txi_refined[[component_name]]
  if (is.matrix(component) && !is.null(rownames(component)) &&
      identical(rownames(component), original_gene_ids)) {
    txi_refined[[component_name]] <- component[keep_gene_ids, , drop = FALSE]
  }
}
if (!identical(rownames(txi_refined$counts), keep_gene_ids)) {
  error_exit("Refined-mask tximport object was not subset correctly")
}
if (nrow(txi_refined$counts) != 101678L) {
  error_exit(paste(
    "Expected 101678 genes after refined masking; found",
    nrow(txi_refined$counts)
  ))
}

dds <- DESeqDataSetFromTximport(
  txi_refined, colData = metadata, design = ~ condition
)
genes_before_low_count_filter <- nrow(dds)
keep_low_count <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep_low_count, ]
genes_after_low_count_filter <- nrow(dds)
if (genes_after_low_count_filter != 24734L) {
  error_exit(paste(
    "Expected 24734 genes after the frozen low-count filter; found",
    genes_after_low_count_filter
  ))
}

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "SG", "S"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id")), drop = FALSE]
res_df$DESeq2_status_host_focused_refined_mask <- ifelse(
  is.na(res_df$padj),
  "NA_padj_not_testable",
  ifelse(
    res_df$padj < 0.01 & res_df$log2FoldChange >= 1,
    "GP_up",
    ifelse(
      res_df$padj < 0.01 & res_df$log2FoldChange <= -1,
      "UP_up",
      "non_DE"
    )
  )
)
write_tsv(res_df, outputs$all_results)

normalized <- as.data.frame(counts(dds, normalized = TRUE))
normalized$gene_id <- rownames(normalized)
normalized <- normalized[, c("gene_id", colnames(txi$counts)), drop = FALSE]
write_tsv(normalized, outputs$normalized_counts)
saveRDS(dds, outputs$dds)

mask_audit <- data.frame(
  gene_id = refined_mask_genes,
  mask_source = ifelse(
    refined_mask_genes %in% expanded_mask_genes,
    "frozen_expanded_nucleotide_mask",
    "additional_coding_denominator_supported_sequence"
  ),
  coding_denominator_review_rules = unname(rule_lookup[refined_mask_genes]),
  present_in_original_tximport = refined_mask_genes %in% original_gene_ids,
  present_in_expanded_mask_results =
    refined_mask_genes %in% expanded_results$gene_id,
  removed_before_low_count_filter = TRUE,
  retained_in_refined_mask_results = refined_mask_genes %in% res_df$gene_id,
  stringsAsFactors = FALSE
)
mask_audit$coding_denominator_review_rules[
  is.na(mask_audit$coding_denominator_review_rules)
] <- ""
if (any(mask_audit$retained_in_refined_mask_results)) {
  error_exit("One or more refined-mask genes were retained in DESeq2 results")
}
write_tsv(mask_audit, outputs$mask_exclusion_audit)

comparison <- merge(
  expanded_results[, needed_expanded_columns],
  res_df[, c(
    "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
    "pvalue", "padj", "DESeq2_status_host_focused_refined_mask"
  )],
  by = "gene_id", all = TRUE, sort = FALSE,
  suffixes = c("_expanded_16112_mask", "_refined_16404_mask")
)
comparison$additional_coding_denominator_masked <-
  comparison$gene_id %in% coding_genes
comparison$comparison_status <- ifelse(
  comparison$additional_coding_denominator_masked,
  "removed_by_refined_coding_denominator_mask",
  ifelse(
    is.na(comparison$DESeq2_status_host_focused_refined_mask),
    "absent_after_refined_mask_low_count_filter",
    ifelse(
      comparison$DESeq2_status_host_focused_expanded_mask ==
        comparison$DESeq2_status_host_focused_refined_mask,
      "status_unchanged",
      "status_changed_after_refit"
    )
  )
)
write_tsv(comparison, outputs$expanded_comparison)

status_counts <- table(res_df$DESeq2_status_host_focused_refined_mask)
count_status <- function(status) {
  if (status %in% names(status_counts)) unname(status_counts[status]) else 0L
}

filtering_summary <- data.frame(
  metric = c(
    "tximport_genes_before_mask",
    "frozen_expanded_mask_genes",
    "additional_coding_denominator_mask_genes",
    "refined_mask_genes_total",
    "additional_candidates_in_expanded_tested_universe",
    "additional_candidates_not_in_expanded_tested_universe",
    "genes_after_refined_mask_before_low_count_filter",
    "genes_after_frozen_low_count_filter",
    "refined_mask_host_focused_total_DEGs",
    "refined_mask_host_focused_GP_up_DEGs",
    "refined_mask_host_focused_UP_up_DEGs",
    "refined_mask_host_focused_non_DE",
    "refined_mask_host_focused_NA_padj_not_testable"
  ),
  value = c(
    length(original_gene_ids),
    length(expanded_mask_genes),
    length(coding_genes),
    length(refined_mask_genes),
    nrow(coding_tested),
    length(coding_genes) - nrow(coding_tested),
    genes_before_low_count_filter,
    genes_after_low_count_filter,
    sum(res_df$DESeq2_status_host_focused_refined_mask %in% c("GP_up", "UP_up")),
    count_status("GP_up"),
    count_status("UP_up"),
    count_status("non_DE"),
    count_status("NA_padj_not_testable")
  ),
  stringsAsFactors = FALSE
)
write_tsv(filtering_summary, outputs$filtering_summary)

parameters <- data.frame(
  parameter = c(
    "analysis_scope",
    "model",
    "contrast",
    "low_count_filter",
    "DEG_threshold",
    "frozen_expanded_mask_genes",
    "additional_coding_denominator_mask_genes",
    "refined_mask_total_genes",
    "mask_policy"
  ),
  value = c(
    "final refined-mask DESeq2 refit; no sequence search or annotation search",
    "tximport gene-level DESeq2 design ~ condition",
    "SG (GP) versus S (UP)",
    "expected count >=10 in at least 3 samples",
    "padj <0.01 and absolute log2FoldChange >=1; no shrinkage",
    "16112",
    "292",
    "16404",
    paste(
      "Freeze existing expanded mask plus coding-denominator candidates meeting",
      "the transcriptome-wide audit rules; no further provenance expansion"
    )
  ),
  stringsAsFactors = FALSE
)
write_tsv(parameters, outputs$parameters)

input_identity <- data.frame(
  input = basename(required_inputs),
  path = required_inputs,
  MD5_before_analysis = input_md5_before,
  stringsAsFactors = FALSE
)
write_tsv(input_identity, outputs$input_identity)
capture.output(sessionInfo(), file = outputs$session)

input_md5_after <- unname(tools::md5sum(required_inputs))
if (!identical(unname(input_md5_before), input_md5_after)) {
  error_exit("One or more inputs changed during DESeq2 analysis")
}

audit_lines <- c(
  paste0("===== ", audit_title, " ====="),
  paste("Discovery root:", discovery_root),
  "",
  "Scope:",
  "  One final host-focused refit using the frozen refined 16,404-gene mask.",
  "  No sequence search, annotation search or candidate-family search was performed.",
  "  All mixed, 3,948-mask and 16,112-mask results were preserved unchanged.",
  "",
  "Frozen model:",
  "  tximport gene-level input; DESeq2 design ~ condition; S (UP) reference; SG (GP) contrast.",
  "  Low-count filter: expected count >=10 in at least 3 samples.",
  "  Reporting threshold: padj <0.01 and absolute log2FoldChange >=1.",
  "  No log2-fold-change shrinkage was applied.",
  "",
  paste("tximport genes before mask:", length(original_gene_ids)),
  paste("Frozen expanded nucleotide mask:", length(expanded_mask_genes)),
  paste("Additional coding-denominator-supported mask:", length(coding_genes)),
  paste("Refined mask total:", length(refined_mask_genes)),
  paste("Additional candidates removed from expanded tested universe:", nrow(coding_tested)),
  paste("Additional candidates already outside expanded tested universe:",
        length(coding_genes) - nrow(coding_tested)),
  paste("Genes after mask before low-count filtering:", genes_before_low_count_filter),
  paste("Genes after frozen low-count filter:", genes_after_low_count_filter),
  "",
  "Refined-mask host-focused DESeq2 statuses:",
  paste("  GP_up:", count_status("GP_up")),
  paste("  UP_up:", count_status("UP_up")),
  paste("  non_DE:", count_status("non_DE")),
  paste("  NA_padj_not_testable:", count_status("NA_padj_not_testable")),
  "",
  "Interpretive boundary:",
  "  The refined mask is sequence-provenance based and does not prove physical RNA source.",
  "  This remains an analysis of the same six libraries, not an independent replication.",
  "  The 16,404-gene mask is frozen; no further provenance-driven mask expansion is planned.",
  "",
  "Input MD5 verification after output generation: PASS",
  "Errors detected: 0",
  paste0(audit_title, ": PASS")
)
writeLines(audit_lines, outputs$audit)

cat(paste(audit_lines, collapse = "\n"), "\n")
for (path in unlist(outputs)) cat("Output:", path, "\n")

