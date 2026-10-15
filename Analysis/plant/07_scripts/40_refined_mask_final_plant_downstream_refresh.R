#!/usr/bin/env Rscript

# Consolidated final downstream refresh for the frozen 16,404-gene plant mask.
#
# One run performs only the required post-DESeq2 work:
#   1. frozen eggNOG annotation join;
#   2. QC, MA and volcano figures using the validated script-27 implementation;
#   3. GO/canonical-KEGG ORA using the validated script-28 implementation;
#   4. current-resolved GO manuscript-table curation;
#   5. frozen plant-defense candidate status overlay.
#
# It does not fit DESeq2, change DEG calls, search sequences or annotations,
# redefine candidate families, or rebuild the optional discovery ranking.

options(stringsAsFactors = FALSE, scipen = 999)

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
scripts_dir <- file.path(discovery_root, "07_scripts")
branch_root <- file.path(discovery_root, "05_host_focused_DESeq2_sensitivity")
results_dir <- file.path(branch_root, "results")
logs_dir <- file.path(branch_root, "logs")
enrichment_dir <- file.path(branch_root, "enrichment")
defense_root <- file.path(discovery_root, "05_plant_defense_audit")
defense_results <- file.path(defense_root, "results")
control_dir <- file.path(discovery_root, "00_project_control")

for (path in c(
  results_dir, logs_dir, enrichment_dir, defense_results, control_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

final_results <- file.path(
  results_dir,
  "plant_gene_DESeq2_GP_vs_UP_host_focused_refined_mask_all.tsv"
)
final_dds <- file.path(
  results_dir, "plant_gene_dds_host_focused_refined_mask.rds"
)
annotation_donor <- file.path(
  plant_root, "08_results_primary",
  "plant_gene_DESeq2_SG_vs_S_all_primary_eggnog_joined.tsv"
)
refined_mask_file <- file.path(
  discovery_root, "04_cross_species_integration", "results",
  "plant_refined_coding_denominator_sequence_review_mask.tsv"
)
frozen_defense_file <- file.path(
  defense_results, "plant_expanded_mask_defense_candidate_gene_summary.tsv"
)
plot_source <- file.path(
  scripts_dir, "27_expanded_mask_host_focused_QC_and_DE_plots.R"
)
enrichment_source <- file.path(
  scripts_dir, "28_expanded_mask_host_focused_GO_KEGG_enrichment.R"
)

joined_output <- file.path(
  results_dir,
  "plant_gene_DESeq2_GP_vs_UP_host_focused_refined_mask_eggnog_joined.tsv"
)
annotation_audit <- file.path(
  logs_dir, "40_refined_mask_host_focused_annotation_join_audit.txt"
)
plot_audit <- file.path(
  logs_dir, "40_refined_mask_host_focused_QC_and_DE_plot_audit.txt"
)
enrichment_audit <- file.path(
  logs_dir, "40_refined_mask_host_focused_GO_KEGG_enrichment_audit.txt"
)
go_curation_summary <- file.path(
  enrichment_dir,
  "host_focused_refined_mask_GO_manuscript_term_curation_summary.tsv"
)
go_curation_audit <- file.path(
  logs_dir, "40_refined_mask_GO_manuscript_term_curation_audit.txt"
)
defense_audit <- file.path(
  defense_root, "40_refined_mask_plant_defense_status_overlay_audit.txt"
)
parameters_output <- file.path(
  control_dir, "40_refined_mask_final_plant_downstream_refresh_parameters.tsv"
)
master_audit <- file.path(
  logs_dir, "40_refined_mask_final_plant_downstream_refresh_audit.txt"
)

master_title <- "REFINED-MASK FINAL PLANT DOWNSTREAM REFRESH"

stop_audit <- function(message) {
  lines <- c(
    paste0("===== ", master_title, " ====="),
    paste0("ERROR: ", message),
    paste0(master_title, ": CHECK REQUIRED")
  )
  writeLines(lines, master_audit)
  stop(message, call. = FALSE)
}

audit_passes <- function(path, pass_line) {
  if (!file.exists(path)) return(FALSE)
  lines <- readLines(path, warn = FALSE)
  length(lines) > 0L && tail(lines, 1L) == pass_line
}

read_tsv_character <- function(path) {
  read.delim(
    path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, colClasses = "character",
    na.strings = NULL, quote = "", comment.char = ""
  )
}

read_tsv_numeric <- function(path) {
  read.delim(
    path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("NA", ""),
    quote = "", comment.char = ""
  )
}

write_tsv <- function(x, path) {
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA"
  )
}

required_inputs <- c(
  final_results, final_dds, annotation_donor, refined_mask_file,
  frozen_defense_file, plot_source, enrichment_source
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop_audit(paste(
    "Missing required input:", paste(missing_inputs, collapse = " | ")
  ))
}

if (audit_passes(master_audit, paste0(master_title, ": PASS"))) {
  stop(
    "Validated final downstream outputs already exist; no files were changed",
    call. = FALSE
  )
}

input_md5_before <- unname(tools::md5sum(required_inputs))

# ---------------------------------------------------------------------------
# 1. Frozen eggNOG annotation join
# ---------------------------------------------------------------------------

annotation_title <- "REFINED-MASK HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN"
annotation_pass <- paste0(annotation_title, ": PASS")

if (!audit_passes(annotation_audit, annotation_pass)) {
  result_char <- read_tsv_character(final_results)
  donor_char <- read_tsv_character(annotation_donor)

  required_result_fields <- c(
    "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
    "pvalue", "padj", "DESeq2_status_host_focused_refined_mask"
  )
  if (!all(required_result_fields %in% names(result_char))) {
    stop_audit(paste(
      "Final DESeq2 result lacks fields:",
      paste(setdiff(required_result_fields, names(result_char)), collapse = ", ")
    ))
  }
  if (nrow(result_char) != 24734L || anyDuplicated(result_char$gene_id)) {
    stop_audit("Final DESeq2 input must contain 24734 unique genes")
  }
  if (nrow(donor_char) != 28189L || anyDuplicated(donor_char$gene_id)) {
    stop_audit("Frozen annotation donor must contain 28189 unique genes")
  }

  donor_index <- match(result_char$gene_id, donor_char$gene_id)
  if (anyNA(donor_index)) {
    stop_audit(paste(
      "Final genes absent from annotation donor:", sum(is.na(donor_index))
    ))
  }

  deseq_fields <- c(
    "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
    "pvalue", "padj", "DESeq2_status",
    "DESeq2_status_host_focused",
    "DESeq2_status_host_focused_expanded_mask",
    "DESeq2_status_host_focused_refined_mask",
    "eggnog_annotation_match"
  )
  annotation_fields <- setdiff(names(donor_char), deseq_fields)
  output_fields <- c(
    names(result_char), "eggnog_annotation_match", annotation_fields
  )
  if (anyDuplicated(output_fields)) {
    stop_audit("Annotation join would create duplicate output columns")
  }

  joined <- result_char
  joined$eggnog_annotation_match <- "matched"
  for (field in annotation_fields) {
    joined[[field]] <- donor_char[[field]][donor_index]
  }
  joined <- joined[, output_fields, drop = FALSE]
  write_tsv(joined, joined_output)

  joined_check <- read_tsv_character(joined_output)
  if (!identical(names(joined_check), output_fields) ||
      nrow(joined_check) != 24734L) {
    stop_audit("Annotation-joined output failed header/row verification")
  }
  for (field in names(result_char)) {
    if (!identical(joined_check[[field]], result_char[[field]])) {
      stop_audit(paste("DESeq2 field changed during annotation join:", field))
    }
  }

  annotation_present <- if ("has_eggnog_annotation" %in% names(joined)) {
    sum(toupper(joined$has_eggnog_annotation) == "TRUE", na.rm = TRUE)
  } else {
    0L
  }
  annotation_lines <- c(
    paste0("===== ", annotation_title, " ====="),
    paste("Refined-mask DESeq2 genes:", nrow(result_char)),
    paste("Frozen eggNOG donor rows:", nrow(donor_char)),
    paste("Genes with annotation-donor match:", nrow(joined)),
    "Genes without annotation-donor match: 0",
    paste("Genes marked with eggNOG annotation:", annotation_present),
    paste("Annotation columns transferred:", length(annotation_fields)),
    "DESeq2 field preservation: PASS",
    "Duplicate output columns: 0",
    "Errors detected: 0",
    annotation_pass
  )
  writeLines(annotation_lines, annotation_audit)
}

# ---------------------------------------------------------------------------
# 2-3. Reuse validated QC/plot and enrichment implementations in memory.
# ---------------------------------------------------------------------------

adapt_expanded_script <- function(path, kind) {
  code <- paste(readLines(path, warn = FALSE), collapse = "\n")
  replacements <- list(
    c("EXPANDED-MASK", "REFINED-MASK"),
    c("Expanded-mask", "Refined-mask"),
    c("expanded-mask", "refined-mask"),
    c("expanded_mask", "refined_mask"),
    c("24931L", "24734L"),
    c("24931", "24734"),
    c("2349L", "2147L"),
    c("2349", "2147"),
    c("2,349", "2,147"),
    c("1386L", "1383L"),
    c("1386", "1383"),
    c("1,386", "1,383"),
    c("21022L", "21032L"),
    c("21022", "21032"),
    c("21,022", "21,032"),
    c("174L", "172L")
  )
  for (pair in replacements) {
    code <- gsub(pair[1], pair[2], code, fixed = TRUE)
  }
  if (kind == "plot") {
    code <- gsub(
      "27_refined_mask", "40_refined_mask", code, fixed = TRUE
    )
  } else if (kind == "enrichment") {
    code <- gsub(
      "28_refined_mask", "40_refined_mask", code, fixed = TRUE
    )
  }
  code
}

plot_pass <- "REFINED-MASK HOST-FOCUSED PLANT QC AND DE PLOT BUILD: PASS"
if (!audit_passes(plot_audit, plot_pass)) {
  plot_code <- adapt_expanded_script(plot_source, "plot")
  eval(parse(text = plot_code), envir = new.env(parent = globalenv()))
}
if (!audit_passes(plot_audit, plot_pass)) {
  stop_audit("Refined-mask QC/DE plot section did not end in PASS")
}

enrichment_pass <- "REFINED-MASK HOST-FOCUSED PLANT GO/KEGG ENRICHMENT: PASS"
if (!audit_passes(enrichment_audit, enrichment_pass)) {
  enrichment_code <- adapt_expanded_script(enrichment_source, "enrichment")
  eval(parse(text = enrichment_code), envir = new.env(parent = globalenv()))
}
if (!audit_passes(enrichment_audit, enrichment_pass)) {
  stop_audit("Refined-mask GO/KEGG enrichment section did not end in PASS")
}

# ---------------------------------------------------------------------------
# 4. Retain only current-resolved significant GO labels for interpretation.
#    Statistical tables and BH values remain unchanged.
# ---------------------------------------------------------------------------

go_curation_title <- "REFINED-MASK GO MANUSCRIPT TERM CURATION"
go_curation_pass <- paste0(go_curation_title, ": PASS")
go_full_paths <- c(
  GP_up = file.path(
    enrichment_dir,
    "host_focused_refined_mask_GO_GP_up_FDR_0.05_named.tsv"
  ),
  UP_up = file.path(
    enrichment_dir,
    "host_focused_refined_mask_GO_UP_up_FDR_0.05_named.tsv"
  )
)
go_resolved_paths <- c(
  GP_up = file.path(
    enrichment_dir,
    "host_focused_refined_mask_GO_GP_up_FDR_0.05_current_resolved.tsv"
  ),
  UP_up = file.path(
    enrichment_dir,
    "host_focused_refined_mask_GO_UP_up_FDR_0.05_current_resolved.tsv"
  )
)
go_top25_paths <- c(
  GP_up = file.path(
    enrichment_dir,
    "host_focused_refined_mask_GO_GP_up_top25_current_resolved.tsv"
  ),
  UP_up = file.path(
    enrichment_dir,
    "host_focused_refined_mask_GO_UP_up_top25_current_resolved.tsv"
  )
)

if (!audit_passes(go_curation_audit, go_curation_pass)) {
  curation_rows <- list()
  for (direction in names(go_full_paths)) {
    full <- read_tsv_numeric(go_full_paths[[direction]])
    if (!all(c("term_id", "GO_term", "ontology", "pvalue", "padj") %in%
             names(full))) {
      stop_audit(paste("GO table lacks required columns:", direction))
    }
    unresolved <- is.na(full$GO_term) |
      full$GO_term == "Unresolved in current GO.db" |
      is.na(full$ontology) | full$ontology == "Unresolved"
    resolved <- full[!unresolved, , drop = FALSE]
    write_tsv(resolved, go_resolved_paths[[direction]])
    write_tsv(head(resolved, 25L), go_top25_paths[[direction]])
    curation_rows[[direction]] <- data.frame(
      direction = direction,
      significant_GO_categories_full = nrow(full),
      significant_unresolved_excluded = sum(unresolved),
      significant_current_resolved_retained = nrow(resolved),
      stringsAsFactors = FALSE
    )
  }
  curation_summary <- do.call(rbind, curation_rows)
  rownames(curation_summary) <- NULL
  write_tsv(curation_summary, go_curation_summary)
  curation_lines <- c(
    paste0("===== ", go_curation_title, " ====="),
    "Scope: post-statistical curation only; no ORA or BH recalculation.",
    "Full statistical enrichment tables were preserved unchanged.",
    "",
    paste(
      "GP_up significant/current-resolved:",
      curation_summary$significant_GO_categories_full[1], "/",
      curation_summary$significant_current_resolved_retained[1]
    ),
    paste(
      "UP_up significant/current-resolved:",
      curation_summary$significant_GO_categories_full[2], "/",
      curation_summary$significant_current_resolved_retained[2]
    ),
    "Original Fisher P-values and BH-adjusted P-values retained unchanged: PASS",
    "Errors detected: 0",
    go_curation_pass
  )
  writeLines(curation_lines, go_curation_audit)
}

# ---------------------------------------------------------------------------
# 5. Frozen defense-candidate membership with final refined-status overlay.
# ---------------------------------------------------------------------------

defense_title <- "REFINED-MASK PLANT-DEFENSE STATUS OVERLAY"
defense_pass <- paste0(defense_title, ": PASS")
defense_outputs <- list(
  gene_summary = file.path(
    defense_results, "plant_refined_mask_defense_candidate_gene_summary.tsv"
  ),
  status_summary = file.path(
    defense_results, "plant_refined_mask_defense_unique_gene_status_summary.tsv"
  ),
  gp_candidates = file.path(
    defense_results, "plant_refined_mask_defense_GP_up_candidates.tsv"
  ),
  up_candidates = file.path(
    defense_results, "plant_refined_mask_defense_UP_up_candidates.tsv"
  ),
  transitions = file.path(
    defense_results,
    "plant_defense_expanded_to_refined_mask_status_transitions.tsv"
  )
)

if (!audit_passes(defense_audit, defense_pass)) {
  frozen <- read_tsv_character(frozen_defense_file)
  refined_res <- read_tsv_character(joined_output)
  refined_mask <- read_tsv_character(refined_mask_file)

  required_frozen <- c(
    "gene_id", "defense_groups", "expanded_mask_host_status"
  )
  required_res <- c(
    "gene_id", "baseMean", "log2FoldChange", "padj",
    "DESeq2_status_host_focused_refined_mask"
  )
  required_mask <- c(
    "gene_id", "additional_coding_denominator_mask_overlap",
    "refined_mask_status"
  )
  if (!all(required_frozen %in% names(frozen)) || nrow(frozen) != 806L ||
      anyDuplicated(frozen$gene_id)) {
    stop_audit("Frozen defense inventory is not the validated 806-gene table")
  }
  if (!all(required_res %in% names(refined_res)) ||
      nrow(refined_res) != 24734L || anyDuplicated(refined_res$gene_id)) {
    stop_audit("Refined result is not the validated 24734-gene table")
  }
  if (!all(required_mask %in% names(refined_mask)) ||
      nrow(refined_mask) != 118082L || anyDuplicated(refined_mask$gene_id)) {
    stop_audit("Refined-mask table is not the validated transcriptome universe")
  }

  refined_excluded <- refined_mask$gene_id[
    refined_mask$refined_mask_status %in% c(
      "review_exclude_frozen_expanded_nucleotide_mask",
      "review_exclude_additional_coding_denominator_supported_sequence"
    )
  ]
  coding_added <- refined_mask$gene_id[
    toupper(refined_mask$additional_coding_denominator_mask_overlap) == "TRUE"
  ]
  if (length(refined_excluded) != 16404L || length(coding_added) != 292L) {
    stop_audit("Refined/coding mask counts differ from 16404/292")
  }

  res_index <- match(frozen$gene_id, refined_res$gene_id)
  in_mask <- frozen$gene_id %in% refined_excluded
  in_coding_addition <- frozen$gene_id %in% coding_added
  final_status <- rep(
    "not_in_refined_mask_host_focused_tested_universe", nrow(frozen)
  )
  final_status[in_mask] <- "excluded_by_refined_sequence_mask"
  final_status[!is.na(res_index)] <- refined_res[[
    "DESeq2_status_host_focused_refined_mask"
  ]][res_index[!is.na(res_index)]]

  updated <- frozen
  updated$refined_mask_host_status <- final_status
  updated$refined_mask_host_baseMean <- ifelse(
    is.na(res_index), "NA", refined_res$baseMean[res_index]
  )
  updated$refined_mask_host_log2FoldChange <- ifelse(
    is.na(res_index), "NA", refined_res$log2FoldChange[res_index]
  )
  updated$refined_mask_host_padj <- ifelse(
    is.na(res_index), "NA", refined_res$padj[res_index]
  )
  updated$refined_sequence_mask_membership <- ifelse(in_mask, "TRUE", "FALSE")
  updated$additional_coding_denominator_mask_membership <- ifelse(
    in_coding_addition, "TRUE", "FALSE"
  )

  final_status_order <- c(
    "GP_up", "UP_up", "non_DE", "NA_padj_not_testable",
    "excluded_by_refined_sequence_mask",
    "not_in_refined_mask_host_focused_tested_universe"
  )
  order_index <- order(
    match(updated$refined_mask_host_status, final_status_order),
    suppressWarnings(as.numeric(updated$refined_mask_host_padj)),
    updated$gene_id,
    na.last = TRUE
  )
  updated <- updated[order_index, , drop = FALSE]
  if (!setequal(updated$gene_id, frozen$gene_id)) {
    stop_audit("Defense candidate membership changed during overlay")
  }
  write_tsv(updated, defense_outputs$gene_summary)
  write_tsv(
    updated[updated$refined_mask_host_status == "GP_up", , drop = FALSE],
    defense_outputs$gp_candidates
  )
  write_tsv(
    updated[updated$refined_mask_host_status == "UP_up", , drop = FALSE],
    defense_outputs$up_candidates
  )

  defense_groups <- c(
    "PR_protein", "jasmonate", "salicylate",
    "peptidase", "peptidase_inhibitor"
  )
  split_groups <- strsplit(updated$defense_groups, ";", fixed = TRUE)
  status_summary_rows <- list()
  row_counter <- 0L
  for (group in defense_groups) {
    group_member <- vapply(
      split_groups,
      function(values) group %in% trimws(values),
      logical(1)
    )
    for (status in final_status_order) {
      row_counter <- row_counter + 1L
      status_summary_rows[[row_counter]] <- data.frame(
        defense_group = group,
        refined_mask_host_status = status,
        unique_candidate_genes = sum(
          group_member & updated$refined_mask_host_status == status
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  write_tsv(
    do.call(rbind, status_summary_rows), defense_outputs$status_summary
  )

  transitions <- as.data.frame(table(
    expanded_16112_mask_host_status = updated$expanded_mask_host_status,
    refined_16404_mask_host_status = updated$refined_mask_host_status
  ), stringsAsFactors = FALSE)
  transitions <- transitions[transitions$Freq > 0L, , drop = FALSE]
  names(transitions)[3] <- "candidate_genes"
  write_tsv(transitions, defense_outputs$transitions)

  status_counts <- table(factor(
    updated$refined_mask_host_status, levels = final_status_order
  ))
  defense_lines <- c(
    paste0("===== ", defense_title, " ====="),
    "Scope: frozen 806-gene candidate membership with final refined-mask status overlay.",
    "No family search, MEROPS search, SignalP/DeepTMHMM run or DESeq2 analysis was performed.",
    "",
    paste("Frozen candidate genes:", nrow(updated)),
    "Candidate membership changed: 0",
    paste("Candidates overlapping complete refined mask:", sum(in_mask)),
    paste("Candidates overlapping 292-gene coding-denominator addition:",
          sum(in_coding_addition)),
    "",
    "Final refined-mask candidate statuses:",
    vapply(
      final_status_order,
      function(status) paste0("  ", status, ": ", status_counts[status]),
      character(1)
    ),
    "",
    "MEROPS absence remains non-negative evidence because the retained screen had restricted scope.",
    "Candidate membership and existing annotation fields preserved: PASS",
    "Errors detected: 0",
    defense_pass
  )
  writeLines(defense_lines, defense_audit)
}

# ---------------------------------------------------------------------------
# Final consolidated accounting
# ---------------------------------------------------------------------------

enrichment_summary <- read_tsv_numeric(file.path(
  enrichment_dir, "host_focused_refined_mask_GO_KEGG_summary.tsv"
))
curation_summary <- read_tsv_numeric(go_curation_summary)
defense_final <- read_tsv_character(defense_outputs$gene_summary)

required_section_passes <- c(
  audit_passes(annotation_audit, annotation_pass),
  audit_passes(plot_audit, plot_pass),
  audit_passes(enrichment_audit, enrichment_pass),
  audit_passes(go_curation_audit, go_curation_pass),
  audit_passes(defense_audit, defense_pass)
)
if (!all(required_section_passes)) {
  stop_audit("One or more consolidated downstream sections did not pass")
}

if (!identical(input_md5_before, unname(tools::md5sum(required_inputs)))) {
  stop_audit("One or more frozen inputs changed during downstream refresh")
}

parameter_table <- data.frame(
  parameter = c(
    "authoritative_tested_genes", "authoritative_GP_up",
    "authoritative_UP_up", "authoritative_total_DEGs",
    "DEG_rule", "enrichment_test", "GO_interpretation_policy",
    "defense_candidate_policy", "discovery_ranking_policy"
  ),
  value = c(
    "24734", "2147", "1383", "3530",
    "padj <0.01 and absolute log2FoldChange >=1; no shrinkage",
    "one-sided Fisher ORA; valid-padj annotated background; category minimum 5; BH <0.05",
    "preserve full statistical tables; interpret only current-resolved labels",
    "freeze 806-gene membership; update expression status only",
    "not rebuilt; optional prioritization branch closed"
  ),
  stringsAsFactors = FALSE
)
write_tsv(parameter_table, parameters_output)

enrichment_lines <- vapply(
  seq_len(nrow(enrichment_summary)),
  function(i) paste0(
    "  ", enrichment_summary$analysis[i], ": ",
    enrichment_summary$significant_categories_BH_lt_0.05[i]
  ),
  character(1)
)
curation_lines <- vapply(
  seq_len(nrow(curation_summary)),
  function(i) paste0(
    "  ", curation_summary$direction[i], ": ",
    curation_summary$significant_current_resolved_retained[i],
    " current-resolved of ",
    curation_summary$significant_GO_categories_full[i],
    " statistically significant"
  ),
  character(1)
)
defense_counts <- table(defense_final$refined_mask_host_status)

master_lines <- c(
  paste0("===== ", master_title, " ====="),
  "Authoritative refined-mask DESeq2 result: 24734 tested genes; 2147 GP-up; 1383 UP-up.",
  "No DESeq2 fitting, DEG reclassification, sequence search, annotation search or candidate-family search was performed.",
  "Optional discovery ranking was not rebuilt.",
  "",
  "Completed sections:",
  "  Frozen eggNOG annotation join: PASS",
  "  QC/PCA/correlation/MA/volcano refresh: PASS",
  "  GO/canonical-KEGG enrichment refresh: PASS",
  "  Current-resolved GO manuscript curation: PASS",
  "  Frozen plant-defense candidate status overlay: PASS",
  "",
  "Significant enrichment categories (BH < 0.05):",
  enrichment_lines,
  "",
  "GO manuscript interpretation sets:",
  curation_lines,
  "",
  "Final defense-candidate statuses:",
  vapply(
    names(defense_counts),
    function(status) paste0("  ", status, ": ", defense_counts[status]),
    character(1)
  ),
  "",
  "Frozen-input MD5 verification after output generation: PASS",
  "Errors detected: 0",
  paste0(master_title, ": PASS")
)
writeLines(master_lines, master_audit)

cat(paste(master_lines, collapse = "\n"), "\n")
cat("Output:", joined_output, "\n")
cat("Output:", parameters_output, "\n")
cat("Output:", master_audit, "\n")
