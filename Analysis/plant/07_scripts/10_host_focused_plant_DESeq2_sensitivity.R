#!/usr/bin/env Rscript

# Host-focused GP-versus-UP DESeq2 sensitivity analysis.
#
# This is a new analysis branch. It excludes only genes in the predeclared
# confirmed multilevel insect-like sequence mask BEFORE size-factor estimation,
# low-count filtering, dispersion fitting and testing. It does not modify any
# original plant result, object, or count file.

options(stringsAsFactors = FALSE)

required_packages <- c("DESeq2", "tximport", "SummarizedExperiment")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(SummarizedExperiment)
})

# Windows RStudio should use the WSL UNC location. Linux execution falls back
# to ~/projects. A custom location can always be supplied through
# QE_PROJECTS_ROOT without editing this script.
projects_root <- Sys.getenv("QE_PROJECTS_ROOT", unset = "")
if (!nzchar(projects_root)) {
  projects_root <- if (.Platform$OS.type == "windows") {
    "//wsl.localhost/Ubuntu/home/dodi/projects"
  } else {
    file.path(Sys.getenv("HOME"), "projects")
  }
}
project_root <- file.path(projects_root, "qe_cross_species_discovery")
plant_root <- file.path(projects_root, "qe_plant_denovo")
output_root <- file.path(project_root, "05_host_focused_DESeq2_sensitivity")
results_dir <- file.path(output_root, "results")
logs_dir <- file.path(output_root, "logs")

input_txi <- file.path(plant_root, "07_deseq2", "plant_gene_txi.rds")
input_metadata <- file.path(plant_root, "01_input", "expression",
                            "plant_deseq2_samples.tsv")
input_mask <- file.path(project_root, "04_cross_species_integration", "results",
                        "plant_confirmed_multilevel_insect_like_sequence_mask.tsv")
input_mixed_results <- file.path(plant_root, "08_results_primary",
                                 "plant_gene_DESeq2_SG_vs_S_all_primary.tsv")

outputs <- list(
  all_results = file.path(results_dir,
                          "plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv"),
  normalized_counts = file.path(results_dir,
                                "plant_gene_normalized_counts_host_focused.tsv"),
  filtering_summary = file.path(results_dir,
                                "plant_host_focused_DESeq2_filtering_summary.tsv"),
  mask_exclusion_audit = file.path(results_dir,
                                   "plant_confirmed_multilevel_mask_exclusion_audit.tsv"),
  mixed_comparison = file.path(results_dir,
                               "plant_host_focused_vs_mixed_primary_DESeq2_comparison.tsv"),
  dds = file.path(results_dir, "plant_gene_dds_host_focused.rds"),
  session = file.path(logs_dir, "plant_host_focused_DESeq2_sessionInfo.txt"),
  audit = file.path(logs_dir, "plant_host_focused_DESeq2_audit.txt")
)

error_exit <- function(message) {
  cat("===== HOST-FOCUSED PLANT DESEQ2 SENSITIVITY ANALYSIS =====\n",
      "ERROR: ", message, "\n",
      "HOST-FOCUSED PLANT DESEQ2 SENSITIVITY ANALYSIS: CHECK REQUIRED\n",
      sep = "", file = stderr())
  quit(status = 1L)
}

if (dir.exists(output_root) &&
    length(list.files(output_root, recursive = TRUE, all.files = FALSE,
                      no.. = TRUE)) > 0L) {
  error_exit(paste(
    "Host-focused output directory already contains files:", output_root,
    "No files were changed. Archive or rename it before rerunning."
  ))
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(input_txi, input_metadata, input_mask, input_mixed_results)) {
  if (!file.exists(path)) error_exit(paste("Missing required input:", path))
}

sha256 <- function(path) {
  unname(tools::md5sum(path)) # file identity check; named explicitly in audit.
}

read_tsv <- function(path, ...) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, ...)
}

write_tsv <- function(object, path, row.names = FALSE) {
  write.table(object, file = path, sep = "\t", quote = FALSE,
              row.names = row.names, col.names = TRUE, na = "")
}

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

metadata <- read_tsv(input_metadata)
if (!all(c("sample_id", "condition") %in% names(metadata))) {
  error_exit("Sample metadata must contain sample_id and condition")
}
if (anyDuplicated(metadata$sample_id) > 0L) error_exit("Duplicate sample_id in metadata")
rownames(metadata) <- metadata$sample_id
if (!identical(colnames(txi$counts), rownames(metadata))) {
  error_exit("tximport sample columns do not exactly match metadata sample order")
}
if (!setequal(unique(metadata$condition), c("S", "SG"))) {
  error_exit("Expected frozen metadata conditions S and SG")
}
metadata$condition <- factor(metadata$condition, levels = c("S", "SG"))

mask <- read_tsv(input_mask)
needed_mask_columns <- c("gene_id", "host_reanalysis_mask_status")
if (!all(needed_mask_columns %in% names(mask))) {
  error_exit("Confirmed mask lacks gene_id or host_reanalysis_mask_status")
}
masked_genes <- unique(mask$gene_id[
  mask$host_reanalysis_mask_status ==
    "review_exclude_confirmed_multilevel_insect_like_sequence"
])
if (length(masked_genes) != 3948L) {
  error_exit(paste("Expected 3948 confirmed multilevel masked genes; found",
                   length(masked_genes)))
}

original_gene_ids <- rownames(txi$counts)
masked_present <- intersect(masked_genes, original_gene_ids)
masked_absent <- setdiff(masked_genes, original_gene_ids)
keep_genes <- setdiff(original_gene_ids, masked_genes)
if (length(masked_present) == 0L) {
  error_exit("No confirmed masked genes overlap the tximport gene universe")
}

# Preserve every gene-indexed tximport component while subsetting before DESeq2.
txi_host <- txi
for (component_name in names(txi_host)) {
  component <- txi_host[[component_name]]
  if (is.matrix(component) && !is.null(rownames(component)) &&
      identical(rownames(component), original_gene_ids)) {
    txi_host[[component_name]] <- component[keep_genes, , drop = FALSE]
  }
}
if (!identical(rownames(txi_host$counts), keep_genes)) {
  error_exit("Masked tximport counts were not subset correctly")
}

dds <- DESeqDataSetFromTximport(txi_host, colData = metadata, design = ~ condition)
pre_filter_genes <- nrow(dds)
keep_low_count <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep_low_count, ]
post_filter_genes <- nrow(dds)
if (post_filter_genes == 0L) error_exit("No genes remained after the frozen low-count filter")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "SG", "S"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id")), drop = FALSE]
res_df$DESeq2_status_host_focused <- ifelse(
  is.na(res_df$padj), "NA_padj_not_testable",
  ifelse(res_df$padj < 0.01 & res_df$log2FoldChange >= 1, "GP_up",
         ifelse(res_df$padj < 0.01 & res_df$log2FoldChange <= -1, "UP_up", "non_DE"))
)
write_tsv(res_df, outputs$all_results)

normalized <- as.data.frame(counts(dds, normalized = TRUE))
normalized$gene_id <- rownames(normalized)
normalized <- normalized[, c("gene_id", colnames(txi$counts)), drop = FALSE]
write_tsv(normalized, outputs$normalized_counts)
saveRDS(dds, outputs$dds)

mask_audit <- data.frame(
  gene_id = masked_genes,
  present_in_original_tximport = masked_genes %in% original_gene_ids,
  removed_before_low_count_filter = masked_genes %in% masked_present,
  retained_for_host_focused_analysis = FALSE,
  stringsAsFactors = FALSE
)
write_tsv(mask_audit, outputs$mask_exclusion_audit)

mixed <- read_tsv(input_mixed_results)
if (!all(c("gene_id", "baseMean", "log2FoldChange", "padj") %in% names(mixed))) {
  error_exit("Mixed baseline DESeq2 table lacks expected result columns")
}
comparison <- merge(
  mixed[, c("gene_id", "baseMean", "log2FoldChange", "padj")],
  res_df[, c("gene_id", "baseMean", "log2FoldChange", "padj",
             "DESeq2_status_host_focused")],
  by = "gene_id", all = TRUE, sort = FALSE, suffixes = c("_mixed", "_host_focused")
)
comparison$DESeq2_status_mixed <- ifelse(
  is.na(comparison$padj_mixed), "NA_padj_not_testable",
  ifelse(comparison$padj_mixed < 0.01 & comparison$log2FoldChange_mixed >= 1, "GP_up",
         ifelse(comparison$padj_mixed < 0.01 & comparison$log2FoldChange_mixed <= -1,
                "UP_up", "non_DE"))
)
comparison$confirmed_multilevel_masked <- comparison$gene_id %in% masked_genes
write_tsv(comparison, outputs$mixed_comparison)

status_counts <- table(res_df$DESeq2_status_host_focused)
summary <- data.frame(
  metric = c(
    "tximport_genes_before_mask",
    "confirmed_multilevel_mask_genes",
    "masked_genes_present_in_tximport",
    "masked_genes_absent_from_tximport",
    "genes_after_mask_before_low_count_filter",
    "genes_after_frozen_low_count_filter",
    "host_focused_total_DEGs",
    "host_focused_GP_up_DEGs",
    "host_focused_UP_up_DEGs",
    "host_focused_non_DE",
    "host_focused_NA_padj_not_testable"
  ),
  value = c(
    length(original_gene_ids),
    length(masked_genes),
    length(masked_present),
    length(masked_absent),
    pre_filter_genes,
    post_filter_genes,
    sum(res_df$DESeq2_status_host_focused %in% c("GP_up", "UP_up")),
    unname(status_counts["GP_up"]),
    unname(status_counts["UP_up"]),
    unname(status_counts["non_DE"]),
    unname(status_counts["NA_padj_not_testable"])
  ),
  stringsAsFactors = FALSE
)
summary$value[is.na(summary$value)] <- 0
write_tsv(summary, outputs$filtering_summary)

input_identity <- data.frame(
  input = c("plant_gene_txi.rds", "plant_deseq2_samples.tsv",
            "confirmed_multilevel_mask.tsv", "mixed_primary_DESeq2.tsv"),
  path = c(input_txi, input_metadata, input_mask, input_mixed_results),
  MD5 = vapply(c(input_txi, input_metadata, input_mask, input_mixed_results),
               sha256, character(1)),
  stringsAsFactors = FALSE
)

audit_lines <- c(
  "===== HOST-FOCUSED PLANT DESEQ2 SENSITIVITY ANALYSIS =====",
  paste("Discovery root:", project_root),
  "",
  "Scope:",
  "  A new sensitivity branch that excludes only the confirmed multilevel insect-like sequence mask",
  "  before DESeq2 normalization, dispersion fitting and testing.",
  "  The original mixed-provenance plant analysis remains unchanged.",
  "",
  "Frozen model:",
  "  tximport gene-level input; DESeq2 design ~ condition; S (UP) reference; SG (GP) contrast.",
  "  Low-count filter: expected count >=10 in at least 3 samples.",
  "  Reporting threshold: padj < 0.01 and absolute log2FoldChange >=1.",
  "  No log2-fold-change shrinkage was applied.",
  "",
  paste("tximport genes before mask:", length(original_gene_ids)),
  paste("Confirmed multilevel masked genes:", length(masked_genes)),
  paste("Masked genes present in tximport:", length(masked_present)),
  paste("Masked genes absent from tximport:", length(masked_absent)),
  paste("Genes after mask before low-count filtering:", pre_filter_genes),
  paste("Genes after frozen low-count filter:", post_filter_genes),
  "",
  "Host-focused DESeq2 statuses:",
  paste("  GP_up:", sum(res_df$DESeq2_status_host_focused == "GP_up")),
  paste("  UP_up:", sum(res_df$DESeq2_status_host_focused == "UP_up")),
  paste("  non_DE:", sum(res_df$DESeq2_status_host_focused == "non_DE")),
  paste("  NA_padj_not_testable:", sum(res_df$DESeq2_status_host_focused == "NA_padj_not_testable")),
  "",
  "Interpretive boundary:",
  "  This is a conservative sequence-provenance sensitivity analysis, not proof of the physical RNA source.",
  "  It does not recover all possible non-plant sequences and must be interpreted alongside the frozen mixed baseline.",
  "  It must not be presented as an independently replicated plant experiment.",
  "",
  "Errors detected: 0",
  "HOST-FOCUSED PLANT DESEQ2 SENSITIVITY ANALYSIS: PASS"
)
writeLines(audit_lines, outputs$audit)
capture.output(sessionInfo(), file = outputs$session)
write_tsv(input_identity, file.path(logs_dir, "host_focused_DESeq2_input_identity.tsv"))

cat(paste(audit_lines, collapse = "\n"), "\n")
for (path in unlist(outputs)) cat("Output:", path, "\n")
