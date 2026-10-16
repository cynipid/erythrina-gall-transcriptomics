# ============================================================================
# Focused insect larva-up MEROPS inhibitor-family shortlist
#
# Purpose:
#   Extract the exact larva-up candidates from the three inhibitor families
#   that are most relevant to the host-interface hypothesis:
#     I04 = serpin family
#     I19 = pacifastin family
#     I01 = Kazal-type inhibitor family
#
# This is an evidence shortlist, NOT a demonstration of inhibitory activity,
# host secretion, or a direct interaction with any plant protease.
#
# Run in RStudio:
# source("//wsl.localhost/Ubuntu/home/dodi/projects/qe_merops_screen/09_insect_larva_up_inhibitor_shortlist.R")
# ============================================================================

options(stringsAsFactors = FALSE)

base_dir <- "//wsl.localhost/Ubuntu/home/dodi/projects/qe_merops_screen"
results_dir <- file.path(base_dir, "results")

if (!dir.exists(results_dir)) {
  stop("Missing MEROPS results directory: ", results_dir)
}

# Use the candidate-level table, not the broader homologue table.  The latter
# is retained as a provenance record but is not the frozen curated shortlist.
candidate_file <- file.path(
  results_dir,
  "insect_MEROPS_classical_secreted_inhibitor_candidates_curated_v2.tsv"
)
if (!file.exists(candidate_file)) {
  stop("Missing curated-v2 insect inhibitor candidate table: ", candidate_file)
}

shortlist_output <- file.path(
  results_dir,
  "insect_larva_up_I01_I04_I19_inhibitor_shortlist_v2.tsv"
)
summary_output <- file.path(
  results_dir,
  "insect_larva_up_I01_I04_I19_inhibitor_shortlist_summary_v2.tsv"
)
audit_output <- file.path(
  results_dir,
  "insect_larva_up_I01_I04_I19_inhibitor_shortlist_audit_v2.txt"
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
  "gene_id", "merops_family", "DESeq2_status", "analysis_class"
)
missing_columns <- setdiff(required_columns, names(candidates))
if (length(missing_columns)) {
  stop(
    "Curated candidate table is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

target_families <- c("I01", "I04", "I19")
family_labels <- c(
  I01 = "Kazal-type inhibitor-family homology",
  I04 = "serpin-family homology",
  I19 = "pacifastin-family homology"
)

shortlist <- candidates[
  candidates$DESeq2_status == "larva_up" &
    candidates$merops_family %in% target_families,
  , drop = FALSE
]

shortlist$interpretation_scope <- unname(family_labels[shortlist$merops_family])
shortlist$interpretation_caution <-
  "MEROPS homology candidate; inhibitory activity, host exposure, and plant target are untested"

# Put the interpretation-critical columns first, while retaining every
# available original evidence field from the curated table.
front_columns <- c(
  "gene_id", "representative_protein", "merops_identifier", "merops_family",
  "merops_reference_organism", "merops_reference_description",
  "reference_evidence_tier", "curation_reason", "percent_identity", "evalue",
  "bitscore", "subject_coverage_percent", "Description", "Preferred_name",
  "PFAMs", "seed_ortholog", "DESeq2_status", "baseMean", "log2FoldChange",
  "padj", "analysis_class", "interpretation_scope", "interpretation_caution"
)
shortlist <- shortlist[, intersect(c(front_columns, names(shortlist)), names(shortlist)),
                       drop = FALSE]
shortlist <- shortlist[order(shortlist$merops_family,
                              -abs(shortlist$log2FoldChange),
                              shortlist$gene_id), , drop = FALSE]

family_order <- data.frame(
  merops_family = target_families,
  family_label = unname(family_labels[target_families]),
  stringsAsFactors = FALSE
)
summary_table <- merge(
  family_order,
  aggregate(
    gene_id ~ merops_family,
    data = shortlist,
    FUN = function(x) length(unique(x))
  ),
  by = "merops_family",
  all.x = TRUE,
  sort = FALSE
)
names(summary_table)[names(summary_table) == "gene_id"] <- "larva_up_candidate_genes"
summary_table$larva_up_candidate_genes[is.na(summary_table$larva_up_candidate_genes)] <- 0L

if (nrow(shortlist)) {
  expression_summary <- aggregate(
    cbind(baseMean, log2FoldChange) ~ merops_family,
    data = shortlist,
    FUN = function(x) c(min = min(x, na.rm = TRUE),
                        median = median(x, na.rm = TRUE),
                        max = max(x, na.rm = TRUE))
  )
  # Convert aggregate's matrix columns to simple, readable scalar columns.
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
  summary_table <- merge(summary_table, expression_summary,
                         by = "merops_family", all.x = TRUE, sort = FALSE)
}
summary_table <- summary_table[match(target_families, summary_table$merops_family), , drop = FALSE]

# The values below were observed in the frozen curated-v2 family table.
expected_counts <- c(I01 = 2L, I04 = 4L, I19 = 3L)
observed_counts <- setNames(summary_table$larva_up_candidate_genes,
                           summary_table$merops_family)

errors <- character()
if (nrow(shortlist) != 9L) {
  errors <- c(errors, "Expected 9 larva-up I01/I04/I19 candidate genes")
}
if (!identical(as.integer(observed_counts[target_families]),
               as.integer(expected_counts[target_families]))) {
  errors <- c(errors, "Larva-up family counts differ from frozen values: I01=2, I04=4, I19=3")
}
if (anyDuplicated(shortlist$gene_id)) {
  errors <- c(errors, "Shortlist contains duplicated gene IDs")
}
if (any(shortlist$analysis_class != "classical_secreted_core")) {
  errors <- c(errors, "Shortlist includes a protein outside classical_secreted_core")
}

write.table(shortlist, shortlist_output, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")
write.table(summary_table, summary_output, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")

report <- c(
  "===== INSECT LARVA-UP INHIBITOR SHORTLIST V2 =====",
  paste("Input curated candidate table:", candidate_file),
  paste("Larva-up I01/I04/I19 genes:", nrow(shortlist)),
  "",
  capture.output(print(summary_table, row.names = FALSE)),
  "",
  paste("Errors detected:", length(errors)),
  if (length(errors)) paste("ERROR:", errors) else "",
  if (length(errors)) {
    "INSECT LARVA-UP INHIBITOR SHORTLIST V2: CHECK REQUIRED"
  } else {
    "INSECT LARVA-UP INHIBITOR SHORTLIST V2: PASS"
  }
)
writeLines(report, audit_output)

cat(paste(report, collapse = "\n"), "\n")
cat("\nOutput:", shortlist_output, "\n")
cat("Output:", summary_output, "\n")
cat("Output:", audit_output, "\n")

if (length(errors)) {
  stop("Insect larva-up inhibitor shortlist failed")
}
