# 12_candidate_gene_inventory_primary.R
# Update candidate-gene DESeq2 status using the official
# primary DEG criterion.
#
# Primary DEG definition:
# padj < 0.01 and |log2FC| >= 1
#
# Candidate families:
# 1. pectin / pectate lyase-like
# 2. cellulase / endoglucanase-like
# 3. chitinase-like

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 1. Read existing full-genome annotation
# ------------------------------------------------------------

annotation <- read.delim(
  "08_results_primary/gene_annotation_table_from_gff.tsv",
  header = TRUE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 2. Read new official DESeq2 results
# ------------------------------------------------------------

all_res <- read.delim(
  "08_results_primary/deseq2_pupa_vs_larva_all_genes.tsv",
  header = TRUE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Build searchable annotation text
# ------------------------------------------------------------

annotation$search_text <- paste(
  annotation$Geneid,
  annotation$Name,
  annotation$gene_symbol,
  annotation$description,
  sep = " "
)


# ------------------------------------------------------------
# 4. Candidate gene-family patterns
# ------------------------------------------------------------

candidate_patterns <- list(
  
  pectin_pectate_lyase =
    "pectin lyase|pectate|polygalacturonase|pectin",
  
  cellulase_endoglucanase =
    "cellulase|endoglucanase|glucanase",
  
  chitinase =
    "chitinase|endochitinase"
)


# ------------------------------------------------------------
# 5. Search full genome annotation
# ------------------------------------------------------------

extract_candidates <- function(
    data,
    family_name,
    pattern
) {
  
  hits <- data[
    grepl(
      pattern,
      data$search_text,
      ignore.case = TRUE
    ),
  ]
  
  if (nrow(hits) > 0) {
    hits$candidate_family <- family_name
  }
  
  return(hits)
}


candidate_list <- lapply(
  names(candidate_patterns),
  function(fam) {
    
    extract_candidates(
      annotation,
      fam,
      candidate_patterns[[fam]]
    )
  }
)


candidate_inventory <- do.call(
  rbind,
  candidate_list
)


# ------------------------------------------------------------
# 6. Keep useful genome-annotation columns
# ------------------------------------------------------------

candidate_inventory <- candidate_inventory[
  ,
  c(
    "candidate_family",
    "Geneid",
    "chromosome",
    "start",
    "end",
    "strand",
    "Name",
    "gene_symbol",
    "description",
    "gene_biotype"
  )
]


# ------------------------------------------------------------
# 7. Merge with new DESeq2 results
# ------------------------------------------------------------

candidate_primary <- merge(
  candidate_inventory,
  all_res[
    ,
    c(
      "Geneid",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    )
  ],
  by = "Geneid",
  all.x = TRUE
)


# ------------------------------------------------------------
# 8. Assign official DESeq2 status
# ------------------------------------------------------------

candidate_primary$DESeq2_status <- "not_in_DESeq2_filtered_results"


# Gene entered DESeq2 analysis but padj is NA

candidate_primary$DESeq2_status[
  !is.na(candidate_primary$baseMean) &
    is.na(candidate_primary$padj)
] <- "tested_but_NA_padj"


# Gene has valid DESeq2 result but does not satisfy
# the primary DEG definition

candidate_primary$DESeq2_status[
  !is.na(candidate_primary$padj)
] <- "tested_not_primary_DEG"


# Pupa-up primary DEG

candidate_primary$DESeq2_status[
  !is.na(candidate_primary$padj) &
    candidate_primary$padj < 0.01 &
    candidate_primary$log2FoldChange >= 1
] <- "pupa_up_primary"


# Larva-up primary DEG

candidate_primary$DESeq2_status[
  !is.na(candidate_primary$padj) &
    candidate_primary$padj < 0.01 &
    candidate_primary$log2FoldChange <= -1
] <- "larva_up_primary"


# ------------------------------------------------------------
# 9. Sort table
# ------------------------------------------------------------

candidate_primary <- candidate_primary[
  order(
    candidate_primary$candidate_family,
    candidate_primary$chromosome,
    candidate_primary$start
  ),
]


# ------------------------------------------------------------
# 10. Save full candidate inventory
# ------------------------------------------------------------

write.table(
  candidate_primary,
  "08_results_primary/candidate_gene_inventory_primary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Create summary
# ------------------------------------------------------------

candidate_summary <- as.data.frame(
  table(
    candidate_primary$candidate_family,
    candidate_primary$DESeq2_status
  )
)

colnames(candidate_summary) <- c(
  "candidate_family",
  "DESeq2_status",
  "count"
)

candidate_summary <- candidate_summary[
  candidate_summary$count > 0,
]


write.table(
  candidate_summary,
  "08_results_primary/candidate_gene_inventory_primary_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Short output
# ------------------------------------------------------------

cat("\nCandidate-gene primary DEG summary\n")

print(candidate_summary)