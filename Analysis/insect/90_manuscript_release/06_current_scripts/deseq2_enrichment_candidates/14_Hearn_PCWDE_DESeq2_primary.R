# 14_Hearn_PCWDE_DESeq2_primary.R
#
# Connect cleaned Hearn-guided PCWDE candidates
# to official insect DESeq2 results.
#
# Primary DEG definition:
# padj < 0.01 and |log2FC| >= 1

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 1. Read cleaned PCWDE candidate inventory
# ------------------------------------------------------------

pcwde <- read.delim(
  "08_results_primary/Hearn_PCWDE_candidate_audit.tsv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Read official DESeq2 results
# ------------------------------------------------------------

deseq <- read.delim(
  "08_results_primary/deseq2_pupa_vs_larva_all_genes.tsv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 3. Merge candidate genes with DESeq2
# ------------------------------------------------------------

pcwde_de <- merge(
  pcwde,
  deseq[
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
# 4. Assign expression category
# ------------------------------------------------------------

pcwde_de$DE_status <- "not_in_DESeq2_filtered_set"

pcwde_de$DE_status[
  !is.na(pcwde_de$baseMean) &
    is.na(pcwde_de$padj)
] <- "tested_but_NA_padj"

pcwde_de$DE_status[
  !is.na(pcwde_de$padj)
] <- "tested_not_primary_DEG"

pcwde_de$DE_status[
  !is.na(pcwde_de$padj) &
    pcwde_de$padj < 0.01 &
    pcwde_de$log2FoldChange >= 1
] <- "pupa_up_primary"

pcwde_de$DE_status[
  !is.na(pcwde_de$padj) &
    pcwde_de$padj < 0.01 &
    pcwde_de$log2FoldChange <= -1
] <- "larva_up_primary"


# ------------------------------------------------------------
# 5. Sort
# ------------------------------------------------------------

pcwde_de <- pcwde_de[
  order(
    pcwde_de$candidate_family,
    pcwde_de$chromosome,
    pcwde_de$start
  ),
]


# ------------------------------------------------------------
# 6. Save full table
# ------------------------------------------------------------

write.table(
  pcwde_de,
  "08_results_primary/Hearn_PCWDE_candidates_with_DESeq2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------

pcwde_summary <- as.data.frame(
  table(
    pcwde_de$candidate_family,
    pcwde_de$DE_status
  )
)

colnames(pcwde_summary) <- c(
  "candidate_family",
  "DE_status",
  "count"
)

pcwde_summary <- pcwde_summary[
  pcwde_summary$count > 0,
]


write.table(
  pcwde_summary,
  "08_results_primary/Hearn_PCWDE_candidates_with_DESeq2_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\nPCWDE candidate expression summary\n")
print(pcwde_summary)