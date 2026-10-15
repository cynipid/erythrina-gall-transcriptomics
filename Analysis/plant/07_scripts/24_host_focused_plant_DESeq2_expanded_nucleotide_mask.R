#!/usr/bin/env Rscript

# Expanded-mask host-focused GP-versus-UP DESeq2 sensitivity analysis.
#
# This script excludes the complete nucleotide-supported Q. erythrinae-like
# review set before normalization, low-count filtering, dispersion fitting and
# testing. It writes explicitly named outputs alongside, but never overwrites,
# the preceding 3,948-gene-mask sensitivity analysis or the original mixed
# plant analysis.

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

# Windows RStudio uses the WSL UNC path. Linux execution falls back to
# ~/projects. QE_PROJECTS_ROOT may be supplied without editing this script.
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
branch_root <- file.path(discovery_root, "05_host_focused_DESeq2_sensitivity")
results_dir <- file.path(branch_root, "results")
logs_dir <- file.path(branch_root, "logs")

input_txi <- file.path(plant_root, "07_deseq2", "plant_gene_txi.rds")
input_metadata <- file.path(
  plant_root, "01_input", "expression", "plant_deseq2_samples.tsv"
)
input_expanded_mask <- file.path(
  discovery_root, "04_cross_species_integration", "results",
  "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv"
)
input_current_host <- file.path(
  results_dir, "plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv"
)

outputs <- list(
  all_results = file.path(
    results_dir,
    "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_all.tsv"
  ),
  normalized_counts = file.path(
    results_dir,
    "plant_gene_normalized_counts_host_focused_expanded_mask.tsv"
  ),
  filtering_summary = file.path(
    results_dir,
    "plant_host_focused_expanded_mask_DESeq2_filtering_summary.tsv"
  ),
  mask_exclusion_audit = file.path(
    results_dir,
    "plant_expanded_nucleotide_mask_exclusion_audit.tsv"
  ),
  current_comparison = file.path(
    results_dir,
    "plant_host_focused_expanded_mask_vs_current_DESeq2_comparison.tsv"
  ),
  dds = file.path(
    results_dir, "plant_gene_dds_host_focused_expanded_mask.rds"
  ),
  input_identity = file.path(
    logs_dir, "24_host_focused_expanded_mask_DESeq2_input_identity.tsv"
  ),
  session = file.path(
    logs_dir, "24_host_focused_expanded_mask_DESeq2_sessionInfo.txt"
  ),
  audit = file.path(
    logs_dir, "24_host_focused_expanded_mask_DESeq2_audit.txt"
  )
)

audit_title <- "EXPANDED-MASK HOST-FOCUSED PLANT DESEQ2 SENSITIVITY ANALYSIS"

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
  input_txi, input_metadata, input_expanded_mask, input_current_host
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
# Legacy metadata labels are used only to reproduce the fitted contrast:
# SG = GP and S = UP.
metadata$condition <- factor(metadata$condition, levels = c("S", "SG"))

mask <- read_tsv(input_expanded_mask)
needed_mask_columns <- c(
  "gene_id", "current_confirmed_mask_overlap",
  "proposed_expanded_mask_status", "nucleotide_provenance_review_class"
)
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

allowed_mask_statuses <- c(
  "review_exclude_existing_confirmed_multilevel_sequence",
  "review_exclude_additional_strong_nucleotide_similarity",
  "retain_no_strong_nucleotide_Qerythrinae_similarity"
)
unexpected_mask_status <- setdiff(
  unique(mask$proposed_expanded_mask_status), allowed_mask_statuses
)
if (length(unexpected_mask_status) > 0L) {
  error_exit(paste(
    "Unexpected proposed mask status:",
    paste(unexpected_mask_status, collapse = ", ")
  ))
}

existing_mask_genes <- unique(mask$gene_id[
  mask$proposed_expanded_mask_status ==
    "review_exclude_existing_confirmed_multilevel_sequence"
])
additional_mask_genes <- unique(mask$gene_id[
  mask$proposed_expanded_mask_status ==
    "review_exclude_additional_strong_nucleotide_similarity"
])
expanded_mask_genes <- union(existing_mask_genes, additional_mask_genes)
if (length(existing_mask_genes) != 3948L) {
  error_exit(paste(
    "Expected 3948 existing confirmed-mask genes; found",
    length(existing_mask_genes)
  ))
}
if (length(additional_mask_genes) != 12164L) {
  error_exit(paste(
    "Expected 12164 additional nucleotide-supported genes; found",
    length(additional_mask_genes)
  ))
}
if (length(expanded_mask_genes) != 16112L) {
  error_exit(paste(
    "Expected 16112 genes in the expanded mask; found",
    length(expanded_mask_genes)
  ))
}
if (!all(mask$current_confirmed_mask_overlap[
      mask$gene_id %in% existing_mask_genes
    ] == "TRUE")) {
  error_exit("Existing confirmed-mask membership is internally inconsistent")
}

current_host <- read_tsv(input_current_host)
needed_current_columns <- c(
  "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
  "pvalue", "padj", "DESeq2_status_host_focused"
)
if (!all(needed_current_columns %in% names(current_host))) {
  error_exit("Current host-focused DESeq2 table lacks required columns")
}
if (nrow(current_host) != 25779L ||
    anyDuplicated(current_host$gene_id) > 0L) {
  error_exit("Current host-focused DESeq2 table must contain 25779 unique genes")
}
expected_current_counts <- c(
  GP_up = 3179L,
  UP_up = 1396L,
  non_DE = 21008L,
  NA_padj_not_testable = 196L
)
observed_current_counts <- table(current_host$DESeq2_status_host_focused)
for (status in names(expected_current_counts)) {
  observed <- if (status %in% names(observed_current_counts)) {
    unname(observed_current_counts[status])
  } else {
    0L
  }
  if (observed != expected_current_counts[status]) {
    error_exit(paste(
      "Current host-focused status mismatch for", status,
      ": observed", observed, "expected", expected_current_counts[status]
    ))
  }
}
if (length(intersect(current_host$gene_id, existing_mask_genes)) > 0L) {
  error_exit("Current host-focused results unexpectedly contain existing-mask genes")
}

additional_current <- current_host[
  current_host$gene_id %in% additional_mask_genes, , drop = FALSE
]
expected_additional_current <- c(
  GP_up = 813L,
  UP_up = 0L,
  non_DE = 14L,
  NA_padj_not_testable = 21L
)
observed_additional_current <- table(
  additional_current$DESeq2_status_host_focused
)
for (status in names(expected_additional_current)) {
  observed <- if (status %in% names(observed_additional_current)) {
    unname(observed_additional_current[status])
  } else {
    0L
  }
  if (observed != expected_additional_current[status]) {
    error_exit(paste(
      "Additional-mask/current-status mismatch for", status,
      ": observed", observed, "expected", expected_additional_current[status]
    ))
  }
}
if (nrow(additional_current) != 848L) {
  error_exit(paste(
    "Expected 848 additional-mask genes in current tested results; found",
    nrow(additional_current)
  ))
}

# Remove the expanded mask before DESeq2 construction and fitting.
keep_gene_ids <- original_gene_ids[!original_gene_ids %in% expanded_mask_genes]
txi_expanded <- txi
for (component_name in names(txi_expanded)) {
  component <- txi_expanded[[component_name]]
  if (is.matrix(component) && !is.null(rownames(component)) &&
      identical(rownames(component), original_gene_ids)) {
    txi_expanded[[component_name]] <- component[
      keep_gene_ids, , drop = FALSE
    ]
  }
}
if (!identical(rownames(txi_expanded$counts), keep_gene_ids)) {
  error_exit("Expanded-mask tximport object was not subset correctly")
}
if (nrow(txi_expanded$counts) != 101970L) {
  error_exit(paste(
    "Expected 101970 genes after expanded masking; found",
    nrow(txi_expanded$counts)
  ))
}

dds <- DESeqDataSetFromTximport(
  txi_expanded, colData = metadata, design = ~ condition
)
genes_before_low_count_filter <- nrow(dds)
keep_low_count <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep_low_count, ]
genes_after_low_count_filter <- nrow(dds)
if (genes_after_low_count_filter != 24931L) {
  error_exit(paste(
    "Expected 24931 genes after the frozen low-count filter; found",
    genes_after_low_count_filter
  ))
}

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "SG", "S"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c(
  "gene_id", setdiff(names(res_df), "gene_id")
), drop = FALSE]
res_df$DESeq2_status_host_focused_expanded_mask <- ifelse(
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
normalized <- normalized[, c(
  "gene_id", colnames(txi$counts)
), drop = FALSE]
write_tsv(normalized, outputs$normalized_counts)
saveRDS(dds, outputs$dds)

mask_lookup <- mask[match(expanded_mask_genes, mask$gene_id), , drop = FALSE]
mask_audit <- data.frame(
  gene_id = expanded_mask_genes,
  mask_source = ifelse(
    expanded_mask_genes %in% existing_mask_genes,
    "existing_confirmed_multilevel_mask",
    "additional_strong_nucleotide_similarity"
  ),
  nucleotide_provenance_review_class =
    mask_lookup$nucleotide_provenance_review_class,
  present_in_original_tximport = expanded_mask_genes %in% original_gene_ids,
  present_in_current_host_focused_results =
    expanded_mask_genes %in% current_host$gene_id,
  removed_before_low_count_filter = TRUE,
  retained_in_expanded_mask_results = expanded_mask_genes %in% res_df$gene_id,
  stringsAsFactors = FALSE
)
if (any(mask_audit$retained_in_expanded_mask_results)) {
  error_exit("One or more expanded-mask genes were retained in DESeq2 results")
}
write_tsv(mask_audit, outputs$mask_exclusion_audit)

comparison <- merge(
  current_host[, c(
    "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
    "pvalue", "padj", "DESeq2_status_host_focused"
  )],
  res_df[, c(
    "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
    "pvalue", "padj", "DESeq2_status_host_focused_expanded_mask"
  )],
  by = "gene_id", all = TRUE, sort = FALSE,
  suffixes = c("_current_3948_mask", "_expanded_16112_mask")
)
comparison$additional_strong_nucleotide_masked <-
  comparison$gene_id %in% additional_mask_genes
comparison$comparison_status <- ifelse(
  comparison$additional_strong_nucleotide_masked,
  "removed_by_expanded_mask",
  ifelse(
    is.na(comparison$DESeq2_status_host_focused_expanded_mask),
    "absent_after_expanded_mask_low_count_filter",
    ifelse(
      comparison$DESeq2_status_host_focused ==
        comparison$DESeq2_status_host_focused_expanded_mask,
      "status_unchanged",
      "status_changed_after_refit"
    )
  )
)
write_tsv(comparison, outputs$current_comparison)

status_counts <- table(res_df$DESeq2_status_host_focused_expanded_mask)
count_status <- function(status) {
  if (status %in% names(status_counts)) unname(status_counts[status]) else 0L
}
filtering_summary <- data.frame(
  metric = c(
    "tximport_genes_before_mask",
    "existing_confirmed_multilevel_mask_genes",
    "additional_strong_nucleotide_mask_genes",
    "expanded_mask_genes_total",
    "expanded_mask_genes_present_in_tximport",
    "genes_after_expanded_mask_before_low_count_filter",
    "genes_after_frozen_low_count_filter",
    "expanded_mask_host_focused_total_DEGs",
    "expanded_mask_host_focused_GP_up_DEGs",
    "expanded_mask_host_focused_UP_up_DEGs",
    "expanded_mask_host_focused_non_DE",
    "expanded_mask_host_focused_NA_padj_not_testable"
  ),
  value = c(
    length(original_gene_ids),
    length(existing_mask_genes),
    length(additional_mask_genes),
    length(expanded_mask_genes),
    sum(expanded_mask_genes %in% original_gene_ids),
    genes_before_low_count_filter,
    genes_after_low_count_filter,
    sum(res_df$DESeq2_status_host_focused_expanded_mask %in%
          c("GP_up", "UP_up")),
    count_status("GP_up"),
    count_status("UP_up"),
    count_status("non_DE"),
    count_status("NA_padj_not_testable")
  ),
  stringsAsFactors = FALSE
)
write_tsv(filtering_summary, outputs$filtering_summary)

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
  "  A final conservative sensitivity refit excluding the complete proposed",
  "  nucleotide-supported Q. erythrinae-like review mask before DESeq2 fitting.",
  "  Earlier mixed and 3,948-gene-mask results were preserved unchanged.",
  "",
  "Frozen model:",
  "  tximport gene-level input; DESeq2 design ~ condition; S (UP) reference; SG (GP) contrast.",
  "  Low-count filter: expected count >=10 in at least 3 samples.",
  "  Reporting threshold: padj <0.01 and absolute log2FoldChange >=1.",
  "  No log2-fold-change shrinkage was applied.",
  "",
  paste("tximport genes before mask:", length(original_gene_ids)),
  paste("Existing confirmed multilevel mask:", length(existing_mask_genes)),
  paste("Additional strong nucleotide mask:", length(additional_mask_genes)),
  paste("Expanded mask total:", length(expanded_mask_genes)),
  paste("Genes after mask before low-count filtering:",
        genes_before_low_count_filter),
  paste("Genes after frozen low-count filter:", genes_after_low_count_filter),
  "",
  "Expanded-mask host-focused DESeq2 statuses:",
  paste("  GP_up:", count_status("GP_up")),
  paste("  UP_up:", count_status("UP_up")),
  paste("  non_DE:", count_status("non_DE")),
  paste("  NA_padj_not_testable:", count_status("NA_padj_not_testable")),
  "",
  "Interpretive boundary:",
  "  The exclusion mask is based on strong sequence provenance evidence.",
  "  It does not prove the physical RNA source or a dissection/laboratory error.",
  "  This remains a sensitivity analysis of the same six plant libraries,",
  "  not an independently replicated experiment.",
  "",
  "Input MD5 verification after output generation: PASS",
  "Errors detected: 0",
  paste0(audit_title, ": PASS")
)
writeLines(audit_lines, outputs$audit)

cat(paste(audit_lines, collapse = "\n"), "\n")
for (path in unlist(outputs)) cat("Output:", path, "\n")
