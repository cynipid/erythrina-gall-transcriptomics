# ============================================================================
# Insect pupa-up MEROPS inhibitor-family audit, curated v2
#
# Purpose:
#   Preserve the complete pupa-up secreted inhibitor-family set and identify
#   the I43, I63, and I08 families that numerically drive the pupa-up signal.
#
# Important interpretation boundary:
#   These are MEROPS homology candidates. The analysis does not establish
#   inhibitory activity, developmental function, plant exposure, or a direct
#   interaction with plant proteins.
#
# Run in RStudio:
# source("//wsl.localhost/Ubuntu/home/dodi/projects/qe_merops_screen/10_insect_pupa_up_inhibitor_audit_v2.R")
# ============================================================================

options(stringsAsFactors = FALSE)

base_dir <- "//wsl.localhost/Ubuntu/home/dodi/projects/qe_merops_screen"
results_dir <- file.path(base_dir, "results")

candidate_file <- file.path(
  results_dir,
  "insect_MEROPS_classical_secreted_inhibitor_candidates_curated_v2.tsv"
)
if (!file.exists(candidate_file)) {
  stop("Missing curated-v2 insect inhibitor candidate table: ", candidate_file)
}

all_output <- file.path(
  results_dir,
  "insect_pupa_up_secreted_inhibitor_candidates_curated_v2.tsv"
)
family_output <- file.path(
  results_dir,
  "insect_pupa_up_secreted_inhibitor_family_summary_curated_v2.tsv"
)
dominant_output <- file.path(
  results_dir,
  "insect_pupa_up_I43_I63_I08_dominant_family_candidates_v2.tsv"
)
audit_output <- file.path(
  results_dir,
  "insect_pupa_up_secreted_inhibitor_audit_v2.txt"
)

candidates <- read.delim(
  candidate_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

required_columns <- c(
  "gene_id", "merops_family", "DESeq2_status", "analysis_class",
  "log2FoldChange", "padj"
)
missing_columns <- setdiff(required_columns, names(candidates))
if (length(missing_columns)) {
  stop("Candidate table is missing required columns: ",
       paste(missing_columns, collapse = ", "))
}

pupa <- candidates[candidates$DESeq2_status == "pupa_up", , drop = FALSE]

# Keep all available evidence, with the manuscript-relevant fields first.
front_columns <- c(
  "gene_id", "representative_protein", "merops_identifier", "merops_family",
  "merops_reference_organism", "merops_reference_description",
  "merops_title", "reference_evidence_tier", "merops_reference_status",
  "evidence_tier", "percent_identity", "evalue", "bitscore",
  "query_span_coverage_percent", "subject_span_coverage_percent",
  "Description", "Preferred_name", "PFAMs", "seed_ortholog",
  "DESeq2_status", "baseMean", "log2FoldChange", "padj", "analysis_class"
)
pupa <- pupa[, intersect(c(front_columns, names(pupa)), names(pupa)), drop = FALSE]
pupa <- pupa[order(pupa$merops_family, -pupa$log2FoldChange, pupa$gene_id),
             , drop = FALSE]

family_summary <- aggregate(
  gene_id ~ merops_family,
  data = pupa,
  FUN = function(x) length(unique(x))
)
names(family_summary)[2] <- "pupa_up_candidate_genes"

if (nrow(pupa)) {
  expression_summary <- aggregate(
    cbind(baseMean, log2FoldChange) ~ merops_family,
    data = pupa,
    FUN = function(x) c(min = min(x, na.rm = TRUE),
                        median = median(x, na.rm = TRUE),
                        max = max(x, na.rm = TRUE))
  )
  expression_summary <- data.frame(
    merops_family = expression_summary$merops_family,
    minimum_baseMean = expression_summary$baseMean[, "min"],
    median_baseMean = expression_summary$baseMean[, "median"],
    maximum_baseMean = expression_summary$baseMean[, "max"],
    minimum_log2FoldChange = expression_summary$log2FoldChange[, "min"],
    median_log2FoldChange = expression_summary$log2FoldChange[, "median"],
    maximum_log2FoldChange = expression_summary$log2FoldChange[, "max"],
    stringsAsFactors = FALSE
  )
  family_summary <- merge(family_summary, expression_summary,
                          by = "merops_family", all.x = TRUE, sort = FALSE)
}
family_summary <- family_summary[order(-family_summary$pupa_up_candidate_genes,
                                       family_summary$merops_family), , drop = FALSE]

# These three families together account for most pupa-up genes. This subset is
# deliberately labelled "dominant", not "host-interface", because the present
# data cannot establish an external plant-facing function.
dominant_families <- c("I43", "I63", "I08")
dominant <- pupa[pupa$merops_family %in% dominant_families, , drop = FALSE]
dominant$interpretation_scope <-
  "Dominant pupa-up inhibitor-family homology; functional role remains unresolved"
dominant <- dominant[order(match(dominant$merops_family, dominant_families),
                           -dominant$log2FoldChange, dominant$gene_id), , drop = FALSE]

observed_dominant <- table(factor(dominant$merops_family,
                                 levels = dominant_families))
expected_dominant <- c(I43 = 28L, I63 = 19L, I08 = 10L)

errors <- character()
if (nrow(pupa) != 68L) {
  errors <- c(errors, "Expected 68 pupa-up classical-secreted inhibitor candidates")
}
if (anyDuplicated(pupa$gene_id)) {
  errors <- c(errors, "Pupa-up table contains duplicated gene IDs")
}
if (any(pupa$analysis_class != "classical_secreted_core")) {
  errors <- c(errors, "Pupa-up table includes a protein outside classical_secreted_core")
}
if (!identical(as.integer(observed_dominant), as.integer(expected_dominant))) {
  errors <- c(errors, "Dominant-family counts differ from frozen values: I43=28, I63=19, I08=10")
}
if (nrow(dominant) != 57L) {
  errors <- c(errors, "Expected 57 genes in the I43/I63/I08 dominant subset")
}

write.table(pupa, all_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(family_summary, family_output, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")
write.table(dominant, dominant_output, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")

report <- c(
  "===== INSECT PUPA-UP SECRETED INHIBITOR AUDIT V2 =====",
  paste("Input curated candidate table:", candidate_file),
  paste("All pupa-up candidates:", nrow(pupa)),
  paste("I43/I63/I08 dominant-family subset:", nrow(dominant)),
  "",
  "Pupa-up candidates by MEROPS family:",
  capture.output(print(family_summary, row.names = FALSE)),
  "",
  paste("Errors detected:", length(errors)),
  if (length(errors)) paste("ERROR:", errors) else "",
  if (length(errors)) {
    "INSECT PUPA-UP SECRETED INHIBITOR AUDIT V2: CHECK REQUIRED"
  } else {
    "INSECT PUPA-UP SECRETED INHIBITOR AUDIT V2: PASS"
  }
)
writeLines(report, audit_output)

cat(paste(report, collapse = "\n"), "\n")
cat("\nOutput:", all_output, "\n")
cat("Output:", family_output, "\n")
cat("Output:", dominant_output, "\n")
cat("Output:", audit_output, "\n")

if (length(errors)) {
  stop("Insect pupa-up secreted inhibitor audit failed")
}
