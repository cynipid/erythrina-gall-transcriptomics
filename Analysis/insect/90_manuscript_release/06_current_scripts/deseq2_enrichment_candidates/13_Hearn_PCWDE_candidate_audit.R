# 13_Hearn_PCWDE_candidate_audit.R
#
# Reference-guided search for plant cell-wall-degrading enzyme
# (PCWDE) candidate families reported by Hearn et al.
#
# This script searches the whole Q. erythrinae genome annotation.
# DESeq2 expression is NOT considered at this stage.

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
# 4. Define Hearn-guided PCWDE signatures
# ------------------------------------------------------------

patterns <- list(
  
  cellulase_endoglucanase =
    paste(
      "cellulase",
      "endoglucanase",
      "glycosyl hydrolase family 5\\b",
      "glycosyl hydrolase family 12\\b",
      "\\bGlyco_hydro_5\\b",
      "\\bGlyco_hydro_12\\b",
      "PF00150",
      "K01179",
      "COG2730",
      "GO:0008810",
      "GO:0030245",
      sep = "|"
    ),
  
  pectin_pectate_lyase =
    paste(
      "pectate lyase",
      "pectin lyase",
      "Pec_lyase_C",
      "PF00544",
      "K01732",
      "COG3866",
      "GO:0047490",
      sep = "|"
    ),
  
  rhamnogalacturonan_lyase =
    paste(
      "rhamnogalacturonan lyase",
      "rhamnogalacturonate lyase",
      "rhamnogalacturonan endolyase",
      "polysaccharide lyase family 4",
      "K18195",
      "ENOG410YE12",
      "GO:0102210",
      sep = "|"
    )
)


# ------------------------------------------------------------
# 5. Search each candidate family
# ------------------------------------------------------------

candidate_list <- list()

for (family in names(patterns)) {
  
  pattern <- patterns[[family]]
  
  gff_ids <- gff$Geneid[
    grepl(
      pattern,
      gff_search,
      ignore.case = TRUE
    )
  ]
  
  eggnog_ids <- eggnog$Geneid[
    grepl(
      pattern,
      eggnog_search,
      ignore.case = TRUE
    )
  ]
  
  all_ids <- unique(
    c(gff_ids, eggnog_ids)
  )
  
  if (length(all_ids) > 0) {
    
    candidate_list[[family]] <- data.frame(
      Geneid = all_ids,
      candidate_family = family,
      GFF_hit = all_ids %in% gff_ids,
      eggNOG_hit = all_ids %in% eggnog_ids,
      stringsAsFactors = FALSE
    )
  }
}


candidates <- do.call(
  rbind,
  candidate_list
)

rownames(candidates) <- NULL


# ------------------------------------------------------------
# 6. Add genome annotation information
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
# 7. Add eggNOG annotation evidence
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
# 8. Sort candidates
# ------------------------------------------------------------

candidates <- candidates[
  order(
    candidates$candidate_family,
    candidates$chromosome,
    candidates$start
  ),
]


# ------------------------------------------------------------
# 9. Save full audit table
# ------------------------------------------------------------

write.table(
  candidates,
  "08_results_primary/Hearn_PCWDE_candidate_audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Create summary
# ------------------------------------------------------------

summary_table <- aggregate(
  Geneid ~ candidate_family,
  data = candidates,
  FUN = function(x) length(unique(x))
)

colnames(summary_table)[2] <- "candidate_genes"


write.table(
  summary_table,
  "08_results_primary/Hearn_PCWDE_candidate_audit_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Short output
# ------------------------------------------------------------

cat("\nHearn-guided PCWDE candidate audit\n")

print(summary_table)

cat("\nEvidence source\n")

print(
  table(
    candidates$candidate_family,
    candidates$GFF_hit,
    candidates$eggNOG_hit
  )
)