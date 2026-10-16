# 15_Hearn_chitinase_candidate_audit.R
#
# Hearn-inspired genome-wide audit of
# Q. erythrinae chitinase / chitinase-like candidates.
#
# This step identifies candidate genes from
# genome + eggNOG annotation only.
#
# DESeq2 expression is NOT considered yet.

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 1. Read genome annotation
# ------------------------------------------------------------

gff <- read.delim(
  "08_results_primary/gene_annotation_table_from_gff.tsv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Read eggNOG annotation
# ------------------------------------------------------------

eggnog <- read.delim(
  "10_functional_annotation/parsed_mappings/Qerythrinae_emapper_annotations_parsed.tsv",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 3. Build searchable annotation text
# ------------------------------------------------------------

gff_search <- paste(
  ifelse(is.na(gff$gene_symbol), "", gff$gene_symbol),
  ifelse(is.na(gff$Name), "", gff$Name),
  ifelse(is.na(gff$description), "", gff$description),
  sep = " | "
)


eggnog_search <- paste(
  ifelse(is.na(eggnog$Description), "", eggnog$Description),
  ifelse(is.na(eggnog$Preferred_name), "", eggnog$Preferred_name),
  ifelse(is.na(eggnog$GOs), "", eggnog$GOs),
  ifelse(is.na(eggnog$KEGG_ko), "", eggnog$KEGG_ko),
  ifelse(is.na(eggnog$CAZy), "", eggnog$CAZy),
  ifelse(is.na(eggnog$PFAMs), "", eggnog$PFAMs),
  ifelse(is.na(eggnog$eggNOG_OGs), "", eggnog$eggNOG_OGs),
  sep = " | "
)


# ------------------------------------------------------------
# 4. Hearn-inspired chitinase signatures
# ------------------------------------------------------------

chitinase_pattern <- paste(
  "\\bchitinase\\b",
  "endochitinase",
  "chitotriosidase",
  "chitinase-like",
  "\\bGlyco_hydro_18\\b",
  "\\bGH18\\b",
  "\\bPF00704\\b",
  "\\bK01183\\b",
  "GO:0004568",
  sep = "|"
)


# ------------------------------------------------------------
# 5. Search GFF annotation
# ------------------------------------------------------------

gff_ids <- gff$Geneid[
  grepl(
    chitinase_pattern,
    gff_search,
    ignore.case = TRUE
  )
]


# ------------------------------------------------------------
# 6. Search eggNOG annotation
# ------------------------------------------------------------

eggnog_ids <- eggnog$Geneid[
  grepl(
    chitinase_pattern,
    eggnog_search,
    ignore.case = TRUE
  )
]


# ------------------------------------------------------------
# 7. Combine candidate Gene IDs
# ------------------------------------------------------------

all_ids <- unique(
  c(
    gff_ids,
    eggnog_ids
  )
)


candidates <- data.frame(
  Geneid = all_ids,
  candidate_family = "GH18_chitinase",
  GFF_hit = all_ids %in% gff_ids,
  eggNOG_hit = all_ids %in% eggnog_ids,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 8. Add genome annotation
# ------------------------------------------------------------

candidates <- merge(
  candidates,
  gff[
    ,
    c(
      "Geneid",
      "chromosome",
      "start",
      "end",
      "strand",
      "gene_symbol",
      "Name",
      "description",
      "gene_biotype"
    )
  ],
  by = "Geneid",
  all.x = TRUE
)


# ------------------------------------------------------------
# 9. Add eggNOG evidence
# ------------------------------------------------------------

eggnog_keep <- eggnog[
  ,
  c(
    "Geneid",
    "seed_ortholog",
    "Description",
    "Preferred_name",
    "GOs",
    "KEGG_ko",
    "CAZy",
    "PFAMs",
    "eggNOG_OGs"
  )
]


colnames(eggnog_keep)[
  colnames(eggnog_keep) == "Description"
] <- "eggNOG_Description"


candidates <- merge(
  candidates,
  eggnog_keep,
  by = "Geneid",
  all.x = TRUE
)


# ------------------------------------------------------------
# 10. Add evidence flags
# ------------------------------------------------------------

candidate_annotation_text <- paste(
  ifelse(is.na(candidates$description), "", candidates$description),
  ifelse(is.na(candidates$eggNOG_Description), "", candidates$eggNOG_Description),
  ifelse(is.na(candidates$Preferred_name), "", candidates$Preferred_name),
  ifelse(is.na(candidates$GOs), "", candidates$GOs),
  ifelse(is.na(candidates$KEGG_ko), "", candidates$KEGG_ko),
  ifelse(is.na(candidates$CAZy), "", candidates$CAZy),
  ifelse(is.na(candidates$PFAMs), "", candidates$PFAMs),
  sep = " | "
)


candidates$GH18_evidence <- grepl(
  "\\bGlyco_hydro_18\\b|\\bGH18\\b|\\bPF00704\\b",
  candidate_annotation_text,
  ignore.case = TRUE
)


candidates$K01183_evidence <- grepl(
  "\\bK01183\\b",
  candidate_annotation_text,
  ignore.case = TRUE
)


candidates$chitinase_name_evidence <- grepl(
  "chitinase|endochitinase|chitotriosidase",
  candidate_annotation_text,
  ignore.case = TRUE
)


# ------------------------------------------------------------
# 11. Sort
# ------------------------------------------------------------

candidates <- candidates[
  order(
    candidates$chromosome,
    candidates$start
  ),
]


# ------------------------------------------------------------
# 12. Save full candidate audit
# ------------------------------------------------------------

write.table(
  candidates,
  "08_results_primary/Hearn_chitinase_candidate_audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Short summary
# ------------------------------------------------------------

cat("\nHearn-inspired chitinase candidate audit\n")

cat(
  "Candidate genes:",
  length(unique(candidates$Geneid)),
  "\n"
)

cat(
  "GFF hits:",
  sum(candidates$GFF_hit),
  "\n"
)

cat(
  "eggNOG hits:",
  sum(candidates$eggNOG_hit),
  "\n"
)

cat(
  "GH18/domain evidence:",
  sum(candidates$GH18_evidence),
  "\n"
)

cat(
  "K01183 evidence:",
  sum(candidates$K01183_evidence),
  "\n"
)

cat(
  "Chitinase-name evidence:",
  sum(candidates$chitinase_name_evidence),
  "\n"
)