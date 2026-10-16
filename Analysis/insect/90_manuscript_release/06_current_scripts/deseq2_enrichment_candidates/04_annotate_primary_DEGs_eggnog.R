# 04_annotate_primary_DEGs_eggnog.R
# Purpose:
# Match the new primary DEGs to the existing eggNOG annotation.
#
# Primary DEG criterion:
# padj < 0.01 and |log2FC| >= 1

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 1. Read new DEG tables
# ------------------------------------------------------------

pupa_up <- read.delim(
  "08_results_primary/DEGs_pupa_up.tsv",
  header = TRUE,
  check.names = FALSE
)

larva_up <- read.delim(
  "08_results_primary/DEGs_larva_up.tsv",
  header = TRUE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 2. Read existing parsed eggNOG annotation
# ------------------------------------------------------------

eggnog <- read.delim(
  "10_functional_annotation/parsed_mappings/Qerythrinae_emapper_annotations_parsed.tsv",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 3. Match DEGs to eggNOG annotation
# ------------------------------------------------------------

pupa_annotated <- merge(
  pupa_up,
  eggnog,
  by = "Geneid",
  all.x = TRUE,
  sort = FALSE
)

larva_annotated <- merge(
  larva_up,
  eggnog,
  by = "Geneid",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 4. Save annotated DEG tables
# ------------------------------------------------------------

write.table(
  pupa_annotated,
  "08_results_primary/DEGs_pupa_up_eggnog_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  larva_annotated,
  "08_results_primary/DEGs_larva_up_eggnog_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Annotation coverage
# ------------------------------------------------------------

pupa_eggnog <- sum(!is.na(pupa_annotated$seed_ortholog))
larva_eggnog <- sum(!is.na(larva_annotated$seed_ortholog))

pupa_GO <- sum(
  !is.na(pupa_annotated$GOs) &
    pupa_annotated$GOs != "-" &
    pupa_annotated$GOs != ""
)

larva_GO <- sum(
  !is.na(larva_annotated$GOs) &
    larva_annotated$GOs != "-" &
    larva_annotated$GOs != ""
)

pupa_KEGG <- sum(
  !is.na(pupa_annotated$KEGG_Pathway) &
    pupa_annotated$KEGG_Pathway != "-" &
    pupa_annotated$KEGG_Pathway != ""
)

larva_KEGG <- sum(
  !is.na(larva_annotated$KEGG_Pathway) &
    larva_annotated$KEGG_Pathway != "-" &
    larva_annotated$KEGG_Pathway != ""
)


# ------------------------------------------------------------
# 6. Print summary
# ------------------------------------------------------------

cat("\n====================================\n")
cat("PRIMARY DEG ANNOTATION SUMMARY\n")
cat("====================================\n")

cat("\nPupa-up total DEGs :", nrow(pupa_annotated), "\n")
cat("Pupa-up eggNOG annotated :", pupa_eggnog, "\n")
cat("Pupa-up with GO :", pupa_GO, "\n")
cat("Pupa-up with KEGG pathway :", pupa_KEGG, "\n")

cat("\nLarva-up total DEGs :", nrow(larva_annotated), "\n")
cat("Larva-up eggNOG annotated :", larva_eggnog, "\n")
cat("Larva-up with GO :", larva_GO, "\n")
cat("Larva-up with KEGG pathway :", larva_KEGG, "\n")