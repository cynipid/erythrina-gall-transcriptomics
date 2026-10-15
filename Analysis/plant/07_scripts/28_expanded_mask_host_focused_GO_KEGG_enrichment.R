#!/usr/bin/env Rscript

# Expanded-mask host-focused E. variegata GO and canonical KEGG ORA.
# This script does not fit DESeq2, alter DEG calls, or overwrite the earlier
# mixed-provenance or 3,948-gene-mask enrichment results.

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(GO.db)
})

projects_root <- if (.Platform$OS.type == "windows") {
  "//wsl.localhost/Ubuntu/home/dodi/projects"
} else {
  file.path(Sys.getenv("HOME"), "projects")
}

discovery_root <- file.path(projects_root, "qe_cross_species_discovery")
branch_root <- file.path(
  discovery_root, "05_host_focused_DESeq2_sensitivity"
)
plant_root <- file.path(projects_root, "qe_plant_denovo")
enrichment_dir <- file.path(branch_root, "enrichment")
logs_dir <- file.path(branch_root, "logs")

dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

result_file <- file.path(
  branch_root, "results",
  "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_eggnog_joined.tsv"
)
go_map_file <- file.path(
  plant_root, "04_enrichment", "primary_FDR001_LFC1", "go",
  "plant_tested_gene_to_GO.tsv"
)
kegg_map_file <- file.path(
  plant_root, "04_enrichment", "primary_FDR001_LFC1", "kegg",
  "plant_tested_gene_to_KEGG_pathway_canonical.tsv"
)
kegg_lookup_file <- file.path(
  plant_root, "04_enrichment", "kegg",
  "KEGG_reference_pathway_lookup.tsv"
)

prefix <- "host_focused_expanded_mask"
outputs <- list(
  go_gp_all = file.path(enrichment_dir, paste0(prefix, "_GO_GP_up_all_terms_named.tsv")),
  go_up_all = file.path(enrichment_dir, paste0(prefix, "_GO_UP_up_all_terms_named.tsv")),
  go_gp_sig = file.path(enrichment_dir, paste0(prefix, "_GO_GP_up_FDR_0.05_named.tsv")),
  go_up_sig = file.path(enrichment_dir, paste0(prefix, "_GO_UP_up_FDR_0.05_named.tsv")),
  go_gp_top25 = file.path(enrichment_dir, paste0(prefix, "_GO_GP_up_top25_named.tsv")),
  go_up_top25 = file.path(enrichment_dir, paste0(prefix, "_GO_UP_up_top25_named.tsv")),
  kegg_gp_all = file.path(enrichment_dir, paste0(prefix, "_KEGG_GP_up_all_pathways_named.tsv")),
  kegg_up_all = file.path(enrichment_dir, paste0(prefix, "_KEGG_UP_up_all_pathways_named.tsv")),
  kegg_gp_sig = file.path(enrichment_dir, paste0(prefix, "_KEGG_GP_up_FDR_0.05_named.tsv")),
  kegg_up_sig = file.path(enrichment_dir, paste0(prefix, "_KEGG_UP_up_FDR_0.05_named.tsv")),
  kegg_gp_top25 = file.path(enrichment_dir, paste0(prefix, "_KEGG_GP_up_top25_named.tsv")),
  kegg_up_top25 = file.path(enrichment_dir, paste0(prefix, "_KEGG_UP_up_top25_named.tsv")),
  foreground = file.path(enrichment_dir, paste0(prefix, "_enrichment_foreground_gene_ids.tsv")),
  summary = file.path(enrichment_dir, paste0(prefix, "_GO_KEGG_summary.tsv")),
  audit_table = file.path(enrichment_dir, paste0(prefix, "_GO_KEGG_audit.tsv")),
  input_identity = file.path(logs_dir, "28_expanded_mask_host_focused_enrichment_input_identity.tsv"),
  session_info = file.path(logs_dir, "28_expanded_mask_host_focused_enrichment_sessionInfo.txt"),
  audit = file.path(logs_dir, "28_expanded_mask_host_focused_GO_KEGG_enrichment_audit.txt")
)

audit_title <- "EXPANDED-MASK HOST-FOCUSED PLANT GO/KEGG ENRICHMENT"

error_exit <- function(message) {
  writeLines(
    c(
      paste0("===== ", audit_title, " ====="),
      paste0("ERROR: ", message),
      paste0(audit_title, ": CHECK REQUIRED")
    ),
    outputs$audit
  )
  stop(message, call. = FALSE)
}

required_files <- c(result_file, go_map_file, kegg_map_file, kegg_lookup_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  error_exit(paste("Missing required input:", paste(missing_files, collapse = " | ")))
}

if (file.exists(outputs$audit)) {
  prior <- readLines(outputs$audit, warn = FALSE)
  if (length(prior) > 0L && tail(prior, 1L) == paste0(audit_title, ": PASS")) {
    stop(
      "Validated expanded-mask enrichment outputs already exist; no files were changed",
      call. = FALSE
    )
  }
}

input_md5_before <- unname(tools::md5sum(required_files))

read_tsv <- function(path) {
  read.delim(
    path, header = TRUE, sep = "\t", check.names = FALSE,
    na.strings = c("NA", ""), quote = "", comment.char = ""
  )
}

write_tsv <- function(x, path) {
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

res <- read_tsv(result_file)
required_result_columns <- c(
  "gene_id", "baseMean", "log2FoldChange", "padj",
  "DESeq2_status_host_focused_expanded_mask"
)
if (!all(required_result_columns %in% names(res))) {
  error_exit(paste(
    "Expanded-mask result columns missing:",
    paste(setdiff(required_result_columns, names(res)), collapse = ", ")
  ))
}
if (nrow(res) != 24931L || anyDuplicated(res$gene_id)) {
  error_exit("Expanded-mask result row count or gene-ID uniqueness differs from the validated input")
}

status_col <- "DESeq2_status_host_focused_expanded_mask"
status_counts <- table(factor(
  res[[status_col]],
  levels = c("GP_up", "UP_up", "non_DE", "NA_padj_not_testable")
))
expected_status_counts <- c(
  GP_up = 2349L, UP_up = 1386L,
  non_DE = 21022L, NA_padj_not_testable = 174L
)
if (!identical(as.integer(status_counts), as.integer(expected_status_counts))) {
  error_exit("Expanded-mask DESeq2 status counts differ from the validated result")
}

padj_numeric <- suppressWarnings(as.numeric(res$padj))
lfc_numeric <- suppressWarnings(as.numeric(res$log2FoldChange))
expected_status <- rep("non_DE", nrow(res))
expected_status[is.na(padj_numeric)] <- "NA_padj_not_testable"
expected_status[!is.na(padj_numeric) & padj_numeric < 0.01 & lfc_numeric >= 1] <- "GP_up"
expected_status[!is.na(padj_numeric) & padj_numeric < 0.01 & lfc_numeric <= -1] <- "UP_up"
if (!identical(as.character(res[[status_col]]), expected_status)) {
  error_exit("Stored expanded-mask statuses do not reproduce the frozen DEG rule")
}

valid_padj_ids <- sort(unique(res$gene_id[!is.na(padj_numeric)]))
gp_ids <- sort(unique(res$gene_id[res[[status_col]] == "GP_up"]))
up_ids <- sort(unique(res$gene_id[res[[status_col]] == "UP_up"]))
if (length(intersect(gp_ids, up_ids)) > 0L) {
  error_exit("GP-up and UP-up foregrounds overlap")
}

clean_mapping <- function(mapping, label) {
  if (!all(c("gene_id", "term_id") %in% names(mapping))) {
    error_exit(paste(label, "mapping lacks gene_id or term_id"))
  }
  mapping <- mapping[
    !is.na(mapping$gene_id) & nzchar(mapping$gene_id) &
      !is.na(mapping$term_id) & nzchar(mapping$term_id) &
      mapping$term_id != "-",
    c("gene_id", "term_id"), drop = FALSE
  ]
  mapping <- unique(mapping)
  mapping[order(mapping$gene_id, mapping$term_id), , drop = FALSE]
}

go_source <- clean_mapping(read_tsv(go_map_file), "GO")
if (
  nrow(go_source) != 684983L ||
    length(unique(go_source$gene_id)) != 11440L ||
    length(unique(go_source$term_id)) != 11253L
) {
  error_exit("Frozen GO mapping source does not match its validated accounting")
}

kegg_source <- clean_mapping(read_tsv(kegg_map_file), "KEGG")
if (
  nrow(kegg_source) != 24683L ||
    length(unique(kegg_source$gene_id)) != 5993L ||
    length(unique(kegg_source$term_id)) != 393L
) {
  error_exit("Frozen KEGG mapping source does not match its validated accounting")
}
if (any(!grepl("^map[0-9]{5}$", kegg_source$term_id))) {
  error_exit("Noncanonical KEGG pathway identifier detected")
}

# Preserve the validated current-pathway correction before deduplication.
kegg_source$term_id[kegg_source$term_id == "map00072"] <- "map00650"
kegg_source$term_id[kegg_source$term_id == "map01130"] <- "map01110"
kegg_source <- unique(kegg_source)
kegg_source <- kegg_source[order(kegg_source$gene_id, kegg_source$term_id), , drop = FALSE]
if (
  nrow(kegg_source) != 24264L ||
    length(unique(kegg_source$gene_id)) != 5993L ||
    length(unique(kegg_source$term_id)) != 391L ||
    any(kegg_source$term_id %in% c("map00072", "map01130"))
) {
  error_exit("Canonical KEGG remapping failed its validated accounting")
}

go_map <- go_source[go_source$gene_id %in% valid_padj_ids, , drop = FALSE]
kegg_map <- kegg_source[kegg_source$gene_id %in% valid_padj_ids, , drop = FALSE]
go_map <- unique(go_map)
kegg_map <- unique(kegg_map)

run_ora <- function(mapping, foreground_ids, min_background_genes = 5L) {
  background_ids <- sort(unique(mapping$gene_id))
  foreground_ids <- sort(intersect(unique(foreground_ids), background_ids))
  if (length(background_ids) == 0L || length(foreground_ids) == 0L) {
    error_exit("An enrichment background or foreground is empty")
  }

  term_sets <- lapply(split(mapping$gene_id, mapping$term_id), unique)
  term_sets <- term_sets[lengths(term_sets) >= min_background_genes]
  if (length(term_sets) == 0L) error_exit("No categories passed the minimum background size")

  rows <- lapply(names(term_sets), function(term_id) {
    term_genes <- term_sets[[term_id]]
    a <- length(intersect(foreground_ids, term_genes))
    b <- length(foreground_ids) - a
    c_value <- length(term_genes) - a
    d <- length(background_ids) - a - b - c_value
    if (any(c(a, b, c_value, d) < 0L)) {
      error_exit(paste("Invalid contingency table for", term_id))
    }
    fisher_result <- fisher.test(
      matrix(c(a, b, c_value, d), nrow = 2L, byrow = TRUE),
      alternative = "greater"
    )
    data.frame(
      term_id = term_id,
      foreground_genes_with_term = a,
      foreground_gene_set_size = length(foreground_ids),
      background_genes_with_term = length(term_genes),
      background_gene_set_size = length(background_ids),
      fold_enrichment = if (a == 0L) 0 else
        (a / length(foreground_ids)) / (length(term_genes) / length(background_ids)),
      odds_ratio = unname(fisher_result$estimate),
      pvalue = fisher_result$p.value
    )
  })

  ans <- do.call(rbind, rows)
  ans$padj <- p.adjust(ans$pvalue, method = "BH")
  ans$significant_BH_0.05 <- !is.na(ans$padj) & ans$padj < 0.05
  ans <- ans[order(ans$padj, ans$pvalue, ans$term_id), , drop = FALSE]
  rownames(ans) <- NULL
  attr(ans, "background_ids") <- background_ids
  attr(ans, "foreground_ids") <- foreground_ids
  ans
}

go_gp <- run_ora(go_map, gp_ids)
go_up <- run_ora(go_map, up_ids)
kegg_gp <- run_ora(kegg_map, gp_ids)
kegg_up <- run_ora(kegg_map, up_ids)

all_go_ids <- sort(unique(c(go_gp$term_id, go_up$term_id)))
go_lookup_raw <- AnnotationDbi::select(
  GO.db::GO.db, keys = all_go_ids, keytype = "GOID",
  columns = c("TERM", "ONTOLOGY")
)
go_lookup_raw <- go_lookup_raw[!duplicated(go_lookup_raw$GOID), , drop = FALSE]
add_go_names <- function(x) {
  idx <- match(x$term_id, go_lookup_raw$GOID)
  x$GO_term <- go_lookup_raw$TERM[idx]
  x$ontology <- go_lookup_raw$ONTOLOGY[idx]
  x$GO_term[is.na(x$GO_term)] <- "Unresolved in current GO.db"
  x$ontology[is.na(x$ontology)] <- "Unresolved"
  x[, c(
    "term_id", "GO_term", "ontology",
    "foreground_genes_with_term", "foreground_gene_set_size",
    "background_genes_with_term", "background_gene_set_size",
    "fold_enrichment", "odds_ratio", "pvalue", "padj",
    "significant_BH_0.05"
  ), drop = FALSE]
}
go_gp <- add_go_names(go_gp)
go_up <- add_go_names(go_up)

kegg_lookup <- read_tsv(kegg_lookup_file)
if (!all(c("KEGG_pathway_ID", "KEGG_pathway_name") %in% names(kegg_lookup))) {
  error_exit("KEGG lookup lacks KEGG_pathway_ID or KEGG_pathway_name")
}
kegg_lookup <- kegg_lookup[!duplicated(kegg_lookup$KEGG_pathway_ID), , drop = FALSE]
add_kegg_names <- function(x) {
  idx <- match(x$term_id, kegg_lookup$KEGG_pathway_ID)
  x$KEGG_pathway_name <- kegg_lookup$KEGG_pathway_name[idx]
  unresolved <- is.na(x$KEGG_pathway_name) | !nzchar(x$KEGG_pathway_name)
  x$KEGG_pathway_name[unresolved] <- "Unresolved in KEGG pathway lookup"
  x[, c(
    "term_id", "KEGG_pathway_name",
    "foreground_genes_with_term", "foreground_gene_set_size",
    "background_genes_with_term", "background_gene_set_size",
    "fold_enrichment", "odds_ratio", "pvalue", "padj",
    "significant_BH_0.05"
  ), drop = FALSE]
}
kegg_gp <- add_kegg_names(kegg_gp)
kegg_up <- add_kegg_names(kegg_up)

tables <- list(GP_GO = go_gp, UP_GO = go_up, GP_KEGG = kegg_gp, UP_KEGG = kegg_up)
for (label in names(tables)) {
  x <- tables[[label]]
  if (anyDuplicated(x$term_id)) error_exit(paste(label, "contains duplicate category IDs"))
  if (any(x$background_genes_with_term < 5L)) {
    error_exit(paste(label, "contains a category below the five-gene background minimum"))
  }
  if (!isTRUE(all.equal(x$padj, p.adjust(x$pvalue, method = "BH"), tolerance = 1e-14))) {
    error_exit(paste(label, "BH-adjusted P-values failed recalculation"))
  }
}
if (any(kegg_gp$KEGG_pathway_name == "Unresolved in KEGG pathway lookup") ||
    any(kegg_up$KEGG_pathway_name == "Unresolved in KEGG pathway lookup")) {
  error_exit("One or more tested KEGG pathway names are unresolved")
}

go_gp_sig <- go_gp[go_gp$significant_BH_0.05, , drop = FALSE]
go_up_sig <- go_up[go_up$significant_BH_0.05, , drop = FALSE]
kegg_gp_sig <- kegg_gp[kegg_gp$significant_BH_0.05, , drop = FALSE]
kegg_up_sig <- kegg_up[kegg_up$significant_BH_0.05, , drop = FALSE]

foreground_table <- rbind(
  data.frame(analysis = "GO", direction = "GP_up", gene_id = sort(intersect(gp_ids, unique(go_map$gene_id)))),
  data.frame(analysis = "GO", direction = "UP_up", gene_id = sort(intersect(up_ids, unique(go_map$gene_id)))),
  data.frame(analysis = "KEGG", direction = "GP_up", gene_id = sort(intersect(gp_ids, unique(kegg_map$gene_id)))),
  data.frame(analysis = "KEGG", direction = "UP_up", gene_id = sort(intersect(up_ids, unique(kegg_map$gene_id))))
)

summary_table <- data.frame(
  analysis = c("GP_GO", "UP_GO", "GP_KEGG", "UP_KEGG"),
  official_direction_genes = c(length(gp_ids), length(up_ids), length(gp_ids), length(up_ids)),
  annotation_eligible_foreground_genes = c(
    go_gp$foreground_gene_set_size[1], go_up$foreground_gene_set_size[1],
    kegg_gp$foreground_gene_set_size[1], kegg_up$foreground_gene_set_size[1]
  ),
  annotation_background_genes = c(
    go_gp$background_gene_set_size[1], go_up$background_gene_set_size[1],
    kegg_gp$background_gene_set_size[1], kegg_up$background_gene_set_size[1]
  ),
  categories_tested_minimum_background_5 = c(
    nrow(go_gp), nrow(go_up), nrow(kegg_gp), nrow(kegg_up)
  ),
  significant_categories_BH_lt_0.05 = c(
    nrow(go_gp_sig), nrow(go_up_sig), nrow(kegg_gp_sig), nrow(kegg_up_sig)
  )
)

audit_table <- data.frame(
  metric = c(
    "expanded_mask_result_rows", "valid_padj_genes", "GP_up_genes", "UP_up_genes",
    "GO_source_pairs", "GO_source_genes", "GO_source_terms",
    "KEGG_source_pairs_before_remap", "KEGG_source_pairs_after_remap",
    "KEGG_source_genes", "KEGG_current_terms",
    "GP_GO_tested", "UP_GO_tested", "GP_KEGG_tested", "UP_KEGG_tested",
    "GP_GO_significant", "UP_GO_significant",
    "GP_KEGG_significant", "UP_KEGG_significant"
  ),
  observed = c(
    nrow(res), length(valid_padj_ids), length(gp_ids), length(up_ids),
    nrow(go_source), length(unique(go_source$gene_id)), length(unique(go_source$term_id)),
    24683L, nrow(kegg_source), length(unique(kegg_source$gene_id)),
    length(unique(kegg_source$term_id)),
    nrow(go_gp), nrow(go_up), nrow(kegg_gp), nrow(kegg_up),
    nrow(go_gp_sig), nrow(go_up_sig), nrow(kegg_gp_sig), nrow(kegg_up_sig)
  )
)

write_tsv(go_gp, outputs$go_gp_all)
write_tsv(go_up, outputs$go_up_all)
write_tsv(go_gp_sig, outputs$go_gp_sig)
write_tsv(go_up_sig, outputs$go_up_sig)
write_tsv(head(go_gp_sig, 25L), outputs$go_gp_top25)
write_tsv(head(go_up_sig, 25L), outputs$go_up_top25)
write_tsv(kegg_gp, outputs$kegg_gp_all)
write_tsv(kegg_up, outputs$kegg_up_all)
write_tsv(kegg_gp_sig, outputs$kegg_gp_sig)
write_tsv(kegg_up_sig, outputs$kegg_up_sig)
write_tsv(head(kegg_gp_sig, 25L), outputs$kegg_gp_top25)
write_tsv(head(kegg_up_sig, 25L), outputs$kegg_up_top25)
write_tsv(foreground_table, outputs$foreground)
write_tsv(summary_table, outputs$summary)
write_tsv(audit_table, outputs$audit_table)

input_identity <- data.frame(
  input_role = c("expanded_mask_DESeq2_eggnog", "GO_mapping", "KEGG_mapping", "KEGG_lookup"),
  path = required_files,
  md5 = input_md5_before
)
write_tsv(input_identity, outputs$input_identity)
writeLines(capture.output(sessionInfo()), outputs$session_info)

input_md5_after <- unname(tools::md5sum(required_files))
if (!identical(input_md5_before, input_md5_after)) {
  error_exit("One or more input files changed during output generation")
}

expected_output_paths <- unlist(outputs[setdiff(names(outputs), "audit")])
missing_outputs <- expected_output_paths[!file.exists(expected_output_paths)]
if (length(missing_outputs) > 0L) {
  error_exit(paste("Expected output missing:", paste(missing_outputs, collapse = " | ")))
}

audit_lines <- c(
  paste0("===== ", audit_title, " ====="),
  paste("Expanded-mask DESeq2 rows:", nrow(res)),
  paste("Genes with valid padj:", length(valid_padj_ids)),
  paste("Official GP-up genes:", length(gp_ids)),
  paste("Official UP-up genes:", length(up_ids)),
  "Background: valid-padj genes with the corresponding GO or canonical KEGG annotation.",
  "Test: one-sided Fisher exact over-representation; minimum five background genes per category.",
  "Multiple testing: BH correction independently for GP-GO, UP-GO, GP-KEGG and UP-KEGG.",
  "Significant enrichment: BH-adjusted P < 0.05.",
  "Canonical KEGG correction: map00072->map00650 and map01130->map01110 before deduplication.",
  "",
  paste("GP GO annotation-eligible foreground:", go_gp$foreground_gene_set_size[1]),
  paste("UP GO annotation-eligible foreground:", go_up$foreground_gene_set_size[1]),
  paste("GO background genes:", go_gp$background_gene_set_size[1]),
  paste("GO categories tested:", nrow(go_gp)),
  paste("Significant GP-up GO categories:", nrow(go_gp_sig)),
  paste("Significant UP-up GO categories:", nrow(go_up_sig)),
  "",
  paste("GP KEGG annotation-eligible foreground:", kegg_gp$foreground_gene_set_size[1]),
  paste("UP KEGG annotation-eligible foreground:", kegg_up$foreground_gene_set_size[1]),
  paste("KEGG background genes:", kegg_gp$background_gene_set_size[1]),
  paste("Canonical KEGG categories tested:", nrow(kegg_gp)),
  paste("Significant GP-up KEGG categories:", nrow(kegg_gp_sig)),
  paste("Significant UP-up KEGG categories:", nrow(kegg_up_sig)),
  "",
  "Earlier enrichment outputs were preserved; no DESeq2 fitting or DEG reclassification was performed.",
  "Input MD5 verification after output generation: PASS",
  "Errors detected: 0",
  paste0(audit_title, ": PASS")
)
writeLines(audit_lines, outputs$audit)

cat(paste(audit_lines, collapse = "\n"), "\n")
cat("\nTop significant GP-up GO categories:\n")
print(head(go_gp_sig[, c("term_id", "GO_term", "ontology", "foreground_genes_with_term", "fold_enrichment", "padj")], 15L), row.names = FALSE)
cat("\nTop significant UP-up GO categories:\n")
print(head(go_up_sig[, c("term_id", "GO_term", "ontology", "foreground_genes_with_term", "fold_enrichment", "padj")], 15L), row.names = FALSE)
cat("\nSignificant GP-up canonical KEGG pathways:\n")
print(kegg_gp_sig[, c("term_id", "KEGG_pathway_name", "foreground_genes_with_term", "fold_enrichment", "padj")], row.names = FALSE)
cat("\nSignificant UP-up canonical KEGG pathways:\n")
print(kegg_up_sig[, c("term_id", "KEGG_pathway_name", "foreground_genes_with_term", "fold_enrichment", "padj")], row.names = FALSE)
for (path in unlist(outputs)) cat("Output:", path, "\n")
