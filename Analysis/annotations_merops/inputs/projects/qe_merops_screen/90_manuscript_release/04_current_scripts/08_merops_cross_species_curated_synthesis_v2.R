# ==========================================================================
# Curated MEROPS cross-species synthesis v2
#
# Run in RStudio after the Python v2 reparser:
# source("//wsl.localhost/Ubuntu/home/dodi/projects/qe_merops_screen/scripts/08_merops_cross_species_curated_synthesis_v2.R")
#
# This script reports three nested insect homology sets and prioritizes plant
# candidates. It does not infer catalytic activity, inhibition, secretion into
# the interacting organism, or physical inhibitor--protease pairs.
# ==========================================================================

options(stringsAsFactors = FALSE)

# 1. Paths -----------------------------------------------------------------

base <- "//wsl.localhost/Ubuntu/home/dodi/projects/qe_merops_screen"
results_dir <- file.path(base, "results")
insect_project <- "//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided"

insect_broad_file <- file.path(
  results_dir,
  "insect_MEROPS_classical_secreted_inhibitor_homologues_curated_v2.tsv"
)
plant_broad_file <- file.path(
  results_dir,
  "plant_GPup_secretome_MEROPS_peptidase_homologues_curated_v2.tsv"
)
insect_integrated_file <- file.path(
  insect_project,
  "11_signalp_tmhmm",
  "results",
  "signalp6_deeptmhmm_deseq2_primary_13455.tsv"
)

summary_output <- file.path(
  results_dir,
  "MEROPS_cross_species_curated_v2_summary.tsv"
)
fisher_output <- file.path(
  results_dir,
  "MEROPS_insect_secreted_inhibitor_direction_sensitivity_v2.tsv"
)
family_output <- file.path(
  results_dir,
  "MEROPS_insect_inhibitor_family_by_DESeq2_curated_v2.tsv"
)
plant_priority_output <- file.path(
  results_dir,
  "plant_GPup_MEROPS_peptidase_priority_curated_v2.tsv"
)
plant_priority_summary_output <- file.path(
  results_dir,
  "plant_GPup_MEROPS_peptidase_priority_summary_v2.tsv"
)
audit_output <- file.path(
  results_dir,
  "MEROPS_cross_species_curated_v2_audit.txt"
)

required_files <- c(
  insect_broad_file,
  plant_broad_file,
  insect_integrated_file
)
if (any(!file.exists(required_files))) {
  missing_files <- required_files[!file.exists(required_files)]
  stop(
    "Missing required v2 input(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

# 2. Helpers ---------------------------------------------------------------

read_tsv <- function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
}

pick_column <- function(data, choices, label) {
  found <- intersect(choices, names(data))
  if (length(found) == 0) {
    stop(label, " not found; expected one of: ", paste(choices, collapse = ", "))
  }
  found[[1]]
}

require_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

is_true <- function(x) {
  as.character(x) %in% c("TRUE", "True", "true", "1")
}

make_direction_test <- function(candidate_ids, core, set_name) {
  candidate_rows <- core[core$gene_id %in% candidate_ids, , drop = FALSE]
  other_rows <- core[!core$gene_id %in% candidate_ids, , drop = FALSE]

  candidate_pupa <- sum(candidate_rows$DESeq2_status == "pupa_up")
  candidate_larva <- sum(candidate_rows$DESeq2_status == "larva_up")
  other_pupa <- sum(other_rows$DESeq2_status == "pupa_up")
  other_larva <- sum(other_rows$DESeq2_status == "larva_up")

  direction_table <- matrix(
    c(candidate_pupa, candidate_larva, other_pupa, other_larva),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      candidate_set = c(set_name, "other_classical_secreted"),
      direction = c("pupa_up", "larva_up")
    )
  )

  if (sum(direction_table[1, ]) == 0 || sum(direction_table[2, ]) == 0) {
    return(
      data.frame(
        candidate_set = set_name,
        candidate_genes = length(candidate_ids),
        candidate_pupa_up = candidate_pupa,
        candidate_larva_up = candidate_larva,
        candidate_pupa_percent_among_DEGs = NA_real_,
        other_secreted_pupa_up = other_pupa,
        other_secreted_larva_up = other_larva,
        other_pupa_percent_among_DEGs = NA_real_,
        odds_ratio = NA_real_,
        confidence_low = NA_real_,
        confidence_high = NA_real_,
        pvalue = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  fisher_result <- fisher.test(direction_table)
  data.frame(
    candidate_set = set_name,
    candidate_genes = length(candidate_ids),
    candidate_pupa_up = candidate_pupa,
    candidate_larva_up = candidate_larva,
    candidate_pupa_percent_among_DEGs =
      100 * candidate_pupa / (candidate_pupa + candidate_larva),
    other_secreted_pupa_up = other_pupa,
    other_secreted_larva_up = other_larva,
    other_pupa_percent_among_DEGs =
      100 * other_pupa / (other_pupa + other_larva),
    odds_ratio = unname(fisher_result$estimate),
    confidence_low = fisher_result$conf.int[[1]],
    confidence_high = fisher_result$conf.int[[2]],
    pvalue = fisher_result$p.value,
    stringsAsFactors = FALSE
  )
}

# 3. Read and normalize inputs --------------------------------------------

insect <- read_tsv(insect_broad_file)
plant <- read_tsv(plant_broad_file)
integrated <- read_tsv(insect_integrated_file)

require_columns(
  insect,
  c(
    "gene_id",
    "merops_family",
    "merops_reference_status",
    "evidence_tier",
    "curated_AB_candidate",
    "DESeq2_status",
    "analysis_class"
  ),
  "Insect v2 broad table"
)
require_columns(
  plant,
  c(
    "gene_id",
    "merops_identifier",
    "merops_family",
    "merops_reference_status",
    "evidence_tier",
    "curated_AB_candidate",
    "seed_taxonomy_code",
    "Description",
    "Preferred_name",
    "PFAMs",
    "DESeq2_status",
    "analysis_class"
  ),
  "Plant v2 broad table"
)

integrated_id_column <- pick_column(
  integrated,
  c("gene_id", "Geneid", "GeneID"),
  "Integrated insect gene ID"
)
integrated_status_column <- pick_column(
  integrated,
  c("DESeq2_status", "deseq2_status", "DE_status", "DESeq2_primary_status"),
  "Integrated insect DESeq2 status"
)
integrated_class_column <- pick_column(
  integrated,
  c("analysis_class", "Analysis_class", "signalp_tmhmm_class"),
  "Integrated insect analysis class"
)

integrated$gene_id <- as.character(integrated[[integrated_id_column]])
integrated$DESeq2_status <- as.character(integrated[[integrated_status_column]])
integrated$analysis_class <- as.character(integrated[[integrated_class_column]])
insect$gene_id <- as.character(insect$gene_id)
plant$gene_id <- as.character(plant$gene_id)

errors <- character()
warnings <- character()

if (nrow(integrated) != 13455) {
  errors <- c(errors, "Integrated insect table does not contain 13,455 rows")
}
if (anyDuplicated(integrated$gene_id) > 0) {
  errors <- c(errors, "Integrated insect gene IDs are not unique")
}
if (anyDuplicated(insect$gene_id) > 0) {
  errors <- c(errors, "Insect v2 broad candidate gene IDs are not unique")
}
if (anyDuplicated(plant$gene_id) > 0) {
  errors <- c(errors, "Plant v2 broad candidate gene IDs are not unique")
}

core <- integrated[
  integrated$analysis_class == "classical_secreted_core",
  ,
  drop = FALSE
]
if (nrow(core) != 2039) {
  errors <- c(errors, "Insect classical-secreted core does not contain 2,039 genes")
}
if (!all(insect$gene_id %in% core$gene_id)) {
  errors <- c(errors, "One or more insect candidates is outside the classical-secreted core")
}

# 4. Nested insect evidence sets and sensitivity tests --------------------

insect_broad_ids <- unique(insect$gene_id)
insect_ab_ids <- unique(insect$gene_id[is_true(insect$curated_AB_candidate)])
insect_a_ids <- unique(
  insect$gene_id[grepl("^A_", as.character(insect$evidence_tier))]
)

if (!all(insect_a_ids %in% insect_ab_ids)) {
  errors <- c(errors, "Insect Tier A set is not a subset of the A+B set")
}
if (!all(insect_ab_ids %in% insect_broad_ids)) {
  errors <- c(errors, "Insect A+B set is not a subset of the broad set")
}

fisher_results <- rbind(
  make_direction_test(
    insect_broad_ids,
    core,
    "broad_inhibitor_family_homology"
  ),
  make_direction_test(
    insect_ab_ids,
    core,
    "curated_A_plus_B_inhibitor_candidates"
  ),
  make_direction_test(
    insect_a_ids,
    core,
    "Tier_A_named_reference_matches"
  )
)
fisher_results$padj_BH_across_three_nested_tests <- p.adjust(
  fisher_results$pvalue,
  method = "BH"
)

# Family-by-expression table for the curated A+B set.
insect_ab <- insect[insect$gene_id %in% insect_ab_ids, , drop = FALSE]
family_levels <- sort(unique(as.character(insect_ab$merops_family)))
status_levels <- c(
  "pupa_up",
  "larva_up",
  "non_DE",
  "NA_padj_not_testable",
  "excluded_low_count"
)
family_grid <- expand.grid(
  merops_family = family_levels,
  DESeq2_status = status_levels,
  stringsAsFactors = FALSE
)
family_counts <- aggregate(
  gene_id ~ merops_family + DESeq2_status,
  data = insect_ab,
  FUN = function(x) length(unique(x))
)
names(family_counts)[names(family_counts) == "gene_id"] <- "candidate_genes"
family_summary <- merge(
  family_grid,
  family_counts,
  by = c("merops_family", "DESeq2_status"),
  all.x = TRUE,
  sort = FALSE
)
family_summary$candidate_genes[is.na(family_summary$candidate_genes)] <- 0L
family_summary <- family_summary[
  order(
    family_summary$merops_family,
    match(family_summary$DESeq2_status, status_levels)
  ),
  ,
  drop = FALSE
]

# 5. Plant evidence prioritization ----------------------------------------

plant_reference_taxa <- c(
  "3847" = "Glycine max",
  "3885" = "Phaseolus vulgaris",
  "3750" = "Malus domestica",
  "3827" = "Cicer arietinum"
)
plant$seed_taxonomy_name <- unname(
  plant_reference_taxa[as.character(plant$seed_taxonomy_code)]
)
plant$plant_reference_seed_supported <- !is.na(plant$seed_taxonomy_name)

annotation_text <- paste(
  plant$Description,
  plant$Preferred_name,
  plant$PFAMs,
  sep = " | "
)
plant$peptidase_annotation_text_consistent <- grepl(
  paste(
    c(
      "peptidase",
      "protease",
      "subtilisin",
      "vacuolar[- ]processing",
      "Peptidase_"
    ),
    collapse = "|"
  ),
  annotation_text,
  ignore.case = TRUE,
  perl = TRUE
)

plant$tier_A_named_reference <- grepl(
  "^A_",
  as.character(plant$evidence_tier)
)
plant$tier_B_unassigned_reference <- grepl(
  "^B_",
  as.character(plant$evidence_tier)
)
plant$tier_C_nonfunctional_reference <- grepl(
  "^C_",
  as.character(plant$evidence_tier)
)

plant$manuscript_priority <- ifelse(
  plant$tier_A_named_reference &
    plant$plant_reference_seed_supported &
    plant$peptidase_annotation_text_consistent,
  "priority_1_named_reference_plus_plant_seed_plus_annotation",
  ifelse(
    plant$tier_B_unassigned_reference &
      plant$plant_reference_seed_supported &
      plant$peptidase_annotation_text_consistent,
    "priority_2_unassigned_reference_plus_plant_seed_plus_annotation",
    ifelse(
      plant$tier_C_nonfunctional_reference,
      "review_explicit_nonpeptidase_reference_homology",
      ifelse(
        !plant$plant_reference_seed_supported,
        "review_unresolved_seed_provenance",
        "review_incomplete_functional_support"
      )
    )
  )
)

plant_priority_summary <- as.data.frame(
  table(plant$manuscript_priority),
  stringsAsFactors = FALSE
)
names(plant_priority_summary) <- c("manuscript_priority", "candidate_genes")
plant_priority_summary <- plant_priority_summary[
  order(-plant_priority_summary$candidate_genes),
  ,
  drop = FALSE
]

# 6. Summary and formal consistency checks --------------------------------

count_status <- function(data, status) {
  sum(as.character(data$DESeq2_status) == status)
}

summary_table <- data.frame(
  section = c(
    rep("insect", 12),
    rep("plant", 12)
  ),
  metric = c(
    "classical_secreted_core",
    "broad_inhibitor_family_homologues",
    "curated_A_plus_B_inhibitor_candidates",
    "Tier_A_named_reference_matches",
    "Tier_C_or_D_homology_only",
    "curated_AB_pupa_up",
    "curated_AB_larva_up",
    "curated_AB_non_DE",
    "curated_AB_NA_padj_not_testable",
    "curated_AB_excluded_low_count",
    "curated_AB_unique_families",
    "nested_direction_tests",
    "GP_up_classical_secretome_screened",
    "broad_peptidase_homologues",
    "curated_A_plus_B_peptidase_candidates",
    "Tier_A_named_reference_matches",
    "Tier_B_unassigned_family_matches",
    "Tier_C_explicit_nonpeptidase_reference_matches",
    "plant_reference_seed_supported",
    "without_plant_reference_seed_support",
    "annotation_text_consistent",
    "priority_1_candidates",
    "priority_2_candidates",
    "review_candidates"
  ),
  count = c(
    nrow(core),
    length(insect_broad_ids),
    length(insect_ab_ids),
    length(insect_a_ids),
    length(setdiff(insect_broad_ids, insect_ab_ids)),
    count_status(insect_ab, "pupa_up"),
    count_status(insect_ab, "larva_up"),
    count_status(insect_ab, "non_DE"),
    count_status(insect_ab, "NA_padj_not_testable"),
    count_status(insect_ab, "excluded_low_count"),
    length(unique(insect_ab$merops_family)),
    nrow(fisher_results),
    549,
    nrow(plant),
    sum(is_true(plant$curated_AB_candidate)),
    sum(plant$tier_A_named_reference),
    sum(plant$tier_B_unassigned_reference),
    sum(plant$tier_C_nonfunctional_reference),
    sum(plant$plant_reference_seed_supported),
    sum(!plant$plant_reference_seed_supported),
    sum(plant$peptidase_annotation_text_consistent),
    sum(grepl("^priority_1", plant$manuscript_priority)),
    sum(grepl("^priority_2", plant$manuscript_priority)),
    sum(grepl("^review_", plant$manuscript_priority))
  ),
  stringsAsFactors = FALSE
)

if (sum(family_summary$candidate_genes) != nrow(insect_ab)) {
  errors <- c(errors, "Insect family-by-status table does not account for the A+B set")
}
if (sum(plant_priority_summary$candidate_genes) != nrow(plant)) {
  errors <- c(errors, "Plant priority summary does not account for the broad set")
}
if (any(plant$tier_C_nonfunctional_reference & is_true(plant$curated_AB_candidate))) {
  errors <- c(errors, "A Tier C plant record was incorrectly included in the A+B set")
}
if (any(insect$evidence_tier == "C_nonfunctional_reference_homology_only" &
        is_true(insect$curated_AB_candidate))) {
  errors <- c(errors, "A Tier C insect record was incorrectly included in the A+B set")
}
if (any(fisher_results$other_secreted_pupa_up < 0) ||
    any(fisher_results$other_secreted_larva_up < 0)) {
  errors <- c(errors, "A Fisher-test background count is negative")
}

# 7. Write outputs ---------------------------------------------------------

write.table(
  summary_table,
  summary_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  fisher_results,
  fisher_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)
write.table(
  family_summary,
  family_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  plant,
  plant_priority_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)
write.table(
  plant_priority_summary,
  plant_priority_summary_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

report <- c(
  "===== MEROPS CROSS-SPECIES CURATED SYNTHESIS V2 =====",
  "",
  "Summary:",
  capture.output(print(summary_table, row.names = FALSE)),
  "",
  "Nested insect direction tests:",
  capture.output(print(fisher_results, row.names = FALSE)),
  "",
  "Plant manuscript-priority summary:",
  capture.output(print(plant_priority_summary, row.names = FALSE)),
  "",
  paste("Warnings:", length(warnings)),
  if (length(warnings) > 0) paste("WARNING:", warnings) else character(),
  paste("Errors detected:", length(errors)),
  if (length(errors) > 0) paste("ERROR:", errors) else character(),
  if (length(errors) == 0) {
    "MEROPS CROSS-SPECIES CURATED SYNTHESIS V2: PASS"
  } else {
    "MEROPS CROSS-SPECIES CURATED SYNTHESIS V2: CHECK REQUIRED"
  }
)
writeLines(report, audit_output)

cat(paste(report, collapse = "\n"), "\n")
cat("\nOutput:", summary_output, "\n")
cat("Output:", fisher_output, "\n")
cat("Output:", family_output, "\n")
cat("Output:", plant_priority_output, "\n")
cat("Output:", plant_priority_summary_output, "\n")
cat("Output:", audit_output, "\n")

if (length(errors) > 0) {
  stop("MEROPS cross-species curated synthesis v2 failed")
}
