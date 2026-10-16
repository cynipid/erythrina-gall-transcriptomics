# ============================================================================
# Final CAP-domain candidate audit for Quadrastichus erythrinae
#
# Run in RStudio with:
# source("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided/07_deseq2/17_Hearn_CAP_domain_audit.R")
#
# This is a conservative, predefined-family audit. It identifies proteins with
# either a PFAM CAP domain or a strong CAP/SCP/antigen-5 annotation, then joins
# the frozen primary DESeq2 and SignalP/DeepTMHMM classifications.
#
# It does not establish effector function or secretion experimentally.
# Official DEG rule: padj < 0.01 and |log2FoldChange| >= 1.
# ============================================================================

options(stringsAsFactors = FALSE)

# 1. Paths
proj <- "//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided"
result_dir <- file.path(proj, "08_results_primary")
signal_dir <- file.path(proj, "11_signalp_tmhmm", "results")
annotation_dir <- file.path(proj, "10_functional_annotation", "eggnog_out")

integrated_file <- file.path(
  signal_dir, "signalp6_deeptmhmm_deseq2_primary_13455.tsv"
)
candidate_output <- file.path(result_dir, "Hearn_CAP_domain_candidate_audit.tsv")
summary_output <- file.path(result_dir, "Hearn_CAP_domain_candidate_audit_summary.tsv")
parameter_output <- file.path(result_dir, "Hearn_CAP_domain_candidate_audit_parameters.tsv")

if (!file.exists(integrated_file)) {
  stop("Missing frozen SignalP/DeepTMHMM-DESeq2 table: ", integrated_file)
}
if (!dir.exists(annotation_dir)) {
  stop("Missing eggNOG output directory: ", annotation_dir)
}

# The representative-protein eggNOG file is selected explicitly by its name.
annotation_files <- list.files(
  annotation_dir, recursive = TRUE, full.names = TRUE,
  pattern = "\\.emapper\\.annotations$"
)
preferred_annotation <- annotation_files[
  grepl("longest.*protein.*per.*gene", basename(annotation_files), ignore.case = TRUE)
]
if (length(preferred_annotation) == 1) {
  annotation_file <- preferred_annotation[[1]]
} else if (length(annotation_files) == 1) {
  # This project uses Qerythrinae.emapper.annotations as the representative
  # protein annotation file; its filename does not contain "longest".
  annotation_file <- annotation_files[[1]]
} else {
  stop(
    "Could not identify exactly one longest-protein eggNOG annotation file. Found:\n",
    paste(annotation_files, collapse = "\n")
  )
}

# 2. Read frozen integrated and representative eggNOG tables
# eggNOG annotation files begin with several ## metadata lines. Read from the
# actual #query header rather than allowing read.delim() to treat metadata as a
# header row.
read_emapper_annotations <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_index <- which(grepl("^#query\\t", lines))[1]
  if (is.na(header_index)) {
    stop("Could not find the #query header in eggNOG annotation file: ", path)
  }
  header <- strsplit(sub("^#", "", lines[[header_index]]), "\t", fixed = TRUE)[[1]]
  data_lines <- lines[seq.int(header_index + 1, length(lines))]
  data_lines <- data_lines[nzchar(data_lines) & !grepl("^#", data_lines)]
  if (length(data_lines) == 0) stop("eggNOG annotation file has no data rows: ", path)
  tab <- read.delim(
    text = data_lines, sep = "\t", header = FALSE, quote = "",
    comment.char = "", fill = TRUE, check.names = FALSE
  )
  if (ncol(tab) < length(header)) {
    stop("eggNOG data rows have fewer columns than the #query header")
  }
  if (ncol(tab) > length(header)) {
    extras <- tab[, seq.int(length(header) + 1, ncol(tab)), drop = FALSE]
    if (any(!is.na(as.matrix(extras)) & nzchar(as.matrix(extras)))) {
      stop("eggNOG data rows have unexpected non-empty columns beyond the #query header")
    }
    tab <- tab[, seq_len(length(header)), drop = FALSE]
  }
  names(tab) <- header
  tab
}

integrated <- read.delim(
  integrated_file, sep = "\t", header = TRUE, quote = "",
  comment.char = "", check.names = FALSE
)
annotation <- read_emapper_annotations(annotation_file)

pick_column <- function(data, choices, label) {
  present <- intersect(choices, names(data))
  if (length(present) == 0) {
    stop(label, " not found. Expected one of: ", paste(choices, collapse = ", "))
  }
  present[[1]]
}

integrated_id_column <- pick_column(
  integrated, c("gene_id", "Geneid", "GeneID"), "Integrated gene-ID column"
)
annotation_id_column <- pick_column(
  annotation, c("#query", "query", "Geneid", "gene_id"), "eggNOG query column"
)

integrated$gene_id <- as.character(integrated[[integrated_id_column]])
annotation$gene_id <- sub("\\s.*$", "", as.character(annotation[[annotation_id_column]]))

# Normalize the two frozen-analysis fields because earlier merge scripts used
# slightly different but equivalent column names.
status_column <- pick_column(
  integrated,
  c("DESeq2_status", "deseq2_status", "DE_status", "DESeq2_primary_status"),
  "Integrated DESeq2-status column"
)
class_column <- pick_column(
  integrated,
  c("analysis_class", "Analysis_class", "signalp_tmhmm_class"),
  "Integrated SignalP/TM-class column"
)
integrated$DESeq2_status <- as.character(integrated[[status_column]])
integrated$analysis_class <- as.character(integrated[[class_column]])

if (nrow(integrated) != 13455) {
  stop("Expected 13,455 integrated protein rows; found ", nrow(integrated))
}
if (anyDuplicated(integrated$gene_id) > 0) {
  stop("Frozen integrated table contains duplicate gene IDs")
}
if (anyDuplicated(annotation$gene_id) > 0) {
  stop("Representative eggNOG table contains duplicate query IDs")
}

# 3. Conservative CAP candidate definition
value_or_blank <- function(data, column) {
  if (column %in% names(data)) as.character(data[[column]]) else rep("", nrow(data))
}

annotation$PFAMs <- value_or_blank(annotation, "PFAMs")
annotation$Description <- value_or_blank(annotation, "Description")
annotation$Preferred_name <- value_or_blank(annotation, "Preferred_name")
annotation$CAZy <- value_or_blank(annotation, "CAZy")
annotation$GOs <- value_or_blank(annotation, "GOs")
annotation$KEGG_ko <- value_or_blank(annotation, "KEGG_ko")
annotation$eggNOG_OGs <- value_or_blank(annotation, "eggNOG_OGs")
annotation$seed_ortholog <- value_or_blank(annotation, "seed_ortholog")

# PFAM CAP is the strongest domain-level criterion. The annotation pattern is
# deliberately narrow to avoid treating generic secretory proteins as CAPs.
annotation$PFAM_CAP <- grepl(
  "(^|[,;[:space:]])CAP($|[,;[:space:]])",
  annotation$PFAMs, perl = TRUE
)
annotation_text <- paste(
  annotation$Description, annotation$Preferred_name, sep = " | "
)
annotation$CAP_annotation <- grepl(
  "CAP (domain|family|superfamily)|cysteine[- ]rich secretory|antigen[- ]?5|venom allergen[- ]?5|SCP/TAPS|pathogenesis[- ]related protein 1",
  annotation_text, ignore.case = TRUE, perl = TRUE
)
annotation$CAP_candidate <- annotation$PFAM_CAP | annotation$CAP_annotation
annotation$CAP_evidence <- ifelse(
  annotation$PFAM_CAP & annotation$CAP_annotation,
  "PFAM_CAP_plus_CAP_annotation",
  ifelse(annotation$PFAM_CAP, "PFAM_CAP", "CAP_annotation_only")
)

cap_annotation <- annotation[annotation$CAP_candidate, , drop = FALSE]
if (anyDuplicated(cap_annotation$gene_id) > 0) {
  stop("CAP annotation subset contains duplicate gene IDs")
}

# 4. Join CAP evidence to the complete frozen integrated protein table
cap_columns <- c(
  "gene_id", "PFAMs", "CAZy", "GOs", "KEGG_ko", "eggNOG_OGs",
  "seed_ortholog", "Description", "Preferred_name", "PFAM_CAP",
  "CAP_annotation", "CAP_evidence"
)
cap_columns <- intersect(cap_columns, names(cap_annotation))
candidates <- merge(
  integrated, cap_annotation[, cap_columns, drop = FALSE],
  by = "gene_id", all = FALSE, sort = FALSE
)

if (nrow(candidates) != nrow(cap_annotation)) {
  stop(
    "Not every CAP candidate was found in the frozen integrated protein table: ",
    nrow(candidates), " of ", nrow(cap_annotation)
  )
}

required_integrated_columns <- c("DESeq2_status", "analysis_class")
missing_integrated_columns <- setdiff(required_integrated_columns, names(candidates))
if (length(missing_integrated_columns) > 0) {
  stop("Merged candidate table is missing: ", paste(missing_integrated_columns, collapse = ", "))
}

front_columns <- c(
  "gene_id", "CAP_evidence", "PFAM_CAP", "CAP_annotation",
  "seed_ortholog", "Description", "Preferred_name", "PFAMs", "CAZy",
  "GOs", "KEGG_ko", "eggNOG_OGs", "baseMean", "log2FoldChange",
  "padj", "DESeq2_status", "analysis_class", "signalp_prediction",
  "deeptmhmm_type", "predicted_TMRs"
)
front_columns <- intersect(front_columns, names(candidates))
candidates <- candidates[, c(front_columns, setdiff(names(candidates), front_columns)), drop = FALSE]

lfc_for_order <- if ("log2FoldChange" %in% names(candidates)) {
  abs(candidates$log2FoldChange)
} else rep(-Inf, nrow(candidates))
candidates <- candidates[order(
  candidates$DESeq2_status != "pupa_up",
  candidates$DESeq2_status != "larva_up",
  -lfc_for_order,
  candidates$gene_id
), , drop = FALSE]

# 5. Family summary
status_levels <- c(
  "pupa_up", "larva_up", "non_DE", "NA_padj_not_testable", "excluded_low_count"
)
class_levels <- c(
  "classical_secreted_core", "signalp_with_downstream_tmr",
  "signalp_tm_discordant_review", "membrane_without_signalp", "no_signalp_no_tmr"
)

summary_table <- data.frame(
  metric = c(
    "CAP_candidate_genes", "PFAM_CAP_candidates", "CAP_annotation_supported_candidates",
    paste0("DESeq2_", status_levels), paste0("analysis_class_", class_levels)
  ),
  count = c(
    nrow(candidates), sum(candidates$PFAM_CAP), sum(candidates$CAP_annotation),
    vapply(status_levels, function(x) sum(candidates$DESeq2_status == x), numeric(1)),
    vapply(class_levels, function(x) sum(candidates$analysis_class == x), numeric(1))
  ),
  stringsAsFactors = FALSE
)

parameters <- data.frame(
  parameter = c(
    "DEG_rule", "candidate_definition", "PFAM_CAP_pattern",
    "CAP_annotation_pattern", "integrated_input", "eggNOG_annotation_input"
  ),
  value = c(
    "padj < 0.01 and |log2FoldChange| >= 1",
    "PFAM CAP domain OR a narrow CAP/SCP/antigen-5 annotation",
    "(^|[,;[:space:]])CAP($|[,;[:space:]])",
    "CAP domain/family/superfamily; cysteine-rich secretory; antigen-5; venom allergen-5; SCP/TAPS; PR-1",
    integrated_file, annotation_file
  ),
  stringsAsFactors = FALSE
)

# 6. Formal audit and output
errors <- character()
if (length(unique(candidates$gene_id)) != nrow(candidates)) {
  errors <- c(errors, "CAP candidate IDs are not unique")
}
if (sum(summary_table$count[grepl("^DESeq2_", summary_table$metric)]) != nrow(candidates)) {
  errors <- c(errors, "DESeq2-status accounting does not equal candidate count")
}
if (sum(summary_table$count[grepl("^analysis_class_", summary_table$metric)]) != nrow(candidates)) {
  errors <- c(errors, "SignalP/TM class accounting does not equal candidate count")
}
if (any(!candidates$PFAM_CAP & !candidates$CAP_annotation)) {
  errors <- c(errors, "One or more reported CAP candidates lacks CAP evidence")
}

write.table(candidates, candidate_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(summary_table, summary_output, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(parameters, parameter_output, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n===== FINAL INSECT CAP-DOMAIN AUDIT =====\n")
cat("Integrated protein-level rows:", nrow(integrated), "\n")
cat("Representative eggNOG annotations:", nrow(annotation), "\n")
cat("CAP candidate genes:", nrow(candidates), "\n\n")
print(summary_table, row.names = FALSE)
cat("\nOutput:", candidate_output, "\n")
cat("Output:", summary_output, "\n")
cat("Output:", parameter_output, "\n")
cat("\nErrors detected:", length(errors), "\n")
if (length(errors) > 0) {
  for (error in errors) cat("ERROR:", error, "\n")
  cat("\nFINAL INSECT CAP-DOMAIN AUDIT: CHECK REQUIRED\n")
  stop("Final insect CAP-domain audit failed")
}
cat("\nFINAL INSECT CAP-DOMAIN AUDIT: PASS\n")
