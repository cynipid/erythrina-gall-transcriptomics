# 16_Hearn_chitinase_DESeq2_primary.R
#
# Connect validated Q. erythrinae chitinase candidates
# to official insect DESeq2 results.
#
# Groups:
# 1. canonical GH18 chitinases
# 2. clean GH19 chitinase class I genes
#
# Primary DEG definition:
# padj < 0.01 and |log2FC| >= 1

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 1. Read chitinase candidate audit
# ------------------------------------------------------------

chit <- read.delim(
  "08_results_primary/Hearn_chitinase_candidate_audit.tsv",
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
# 3. Define canonical GH18 chitinases
# ------------------------------------------------------------

GH18 <- chit[
  grepl(
    "Glyco_hydro_18",
    chit$PFAMs,
    fixed = TRUE
  ) &
    grepl(
      "K01183",
      chit$KEGG_ko,
      fixed = TRUE
    ),
]

GH18$chitinase_group <- "GH18_canonical"


# ------------------------------------------------------------
# 4. Define clean GH19 chitinase class I genes
# ------------------------------------------------------------

GH19 <- chit[
  grepl(
    "Glyco_hydro_19",
    chit$PFAMs,
    fixed = TRUE
  ) &
    grepl(
      "Chitinase class I",
      chit$eggNOG_Description,
      fixed = TRUE
    ),
]

GH19$chitinase_group <- "GH19_class_I"


# ------------------------------------------------------------
# 5. Combine groups
# ------------------------------------------------------------

chitinase_clean <- rbind(
  GH18,
  GH19
)


# ------------------------------------------------------------
# 6. Merge with DESeq2
# ------------------------------------------------------------

chitinase_de <- merge(
  chitinase_clean,
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
# 7. Assign DESeq2 status
# ------------------------------------------------------------

chitinase_de$DE_status <- "not_in_DESeq2_filtered_set"

chitinase_de$DE_status[
  !is.na(chitinase_de$baseMean) &
    is.na(chitinase_de$padj)
] <- "tested_but_NA_padj"

chitinase_de$DE_status[
  !is.na(chitinase_de$padj)
] <- "tested_not_primary_DEG"

chitinase_de$DE_status[
  !is.na(chitinase_de$padj) &
    chitinase_de$padj < 0.01 &
    chitinase_de$log2FoldChange >= 1
] <- "pupa_up_primary"

chitinase_de$DE_status[
  !is.na(chitinase_de$padj) &
    chitinase_de$padj < 0.01 &
    chitinase_de$log2FoldChange <= -1
] <- "larva_up_primary"


# ------------------------------------------------------------
# 8. Sort
# ------------------------------------------------------------

chitinase_de <- chitinase_de[
  order(
    chitinase_de$chitinase_group,
    chitinase_de$chromosome,
    chitinase_de$start
  ),
]


# ------------------------------------------------------------
# 9. Save full result table
# ------------------------------------------------------------

write.table(
  chitinase_de,
  "08_results_primary/Hearn_chitinase_candidates_with_DESeq2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Create summary
# ------------------------------------------------------------

chitinase_summary <- as.data.frame(
  table(
    chitinase_de$chitinase_group,
    chitinase_de$DE_status
  )
)

colnames(chitinase_summary) <- c(
  "chitinase_group",
  "DE_status",
  "count"
)

chitinase_summary <- chitinase_summary[
  chitinase_summary$count > 0,
]


write.table(
  chitinase_summary,
  "08_results_primary/Hearn_chitinase_candidates_with_DESeq2_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Short output
# ------------------------------------------------------------

cat("\nChitinase candidate expression summary\n")

print(chitinase_summary)