# 02_run_deseq2_primary_FDR001_LFC1.R
# Purpose:
# Run DESeq2 differential expression analysis for Q. erythrinae
# Comparison: pupa vs larva
#
# Primary DEG definition:
# padj < 0.01 and |log2FoldChange| >= 1
#
# Interpretation:
# positive log2FC = higher expression in pupa
# negative log2FC = higher expression in larva


# ------------------------------------------------------------
# 1. Set project folder
# ------------------------------------------------------------

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 2. Load DESeq2
# ------------------------------------------------------------

library(DESeq2)


# ------------------------------------------------------------
# 3. Read count matrix
# ------------------------------------------------------------

count_matrix <- read.delim(
  "06_counts/gene_counts_matrix.tsv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Read sample metadata
# ------------------------------------------------------------

metadata <- read.delim(
  "06_counts/sample_metadata.tsv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


# ------------------------------------------------------------
# 5. Check sample order
# ------------------------------------------------------------

print(colnames(count_matrix))
print(rownames(metadata))

stopifnot(all(colnames(count_matrix) == rownames(metadata)))


# ------------------------------------------------------------
# 6. Define stage order
# ------------------------------------------------------------

# Larva = reference condition
# Therefore:
# positive log2FoldChange = pupa higher
# negative log2FoldChange = larva higher

metadata$stage <- factor(
  metadata$stage,
  levels = c("larva", "pupa")
)


# ------------------------------------------------------------
# 7. Create DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ stage
)


# ------------------------------------------------------------
# 8. Low-count filtering
# ------------------------------------------------------------

# Keep genes with total count >= 10 across all 6 samples

dds <- dds[rowSums(counts(dds)) >= 10, ]

cat("\nGenes retained after low-count filtering:",
    nrow(dds), "\n")


# ------------------------------------------------------------
# 9. Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(dds)


# ------------------------------------------------------------
# 10. Extract pupa vs larva results
# ------------------------------------------------------------

res <- results(
  dds,
  contrast = c("stage", "pupa", "larva"),
  alpha = 0.01
)


# ------------------------------------------------------------
# 11. Order results by adjusted P-value
# ------------------------------------------------------------

res_ordered <- res[order(res$padj), ]

res_df <- as.data.frame(res_ordered)

res_df$Geneid <- rownames(res_df)

res_df <- res_df[
  , c("Geneid", setdiff(colnames(res_df), "Geneid"))
]


# ------------------------------------------------------------
# 12. Create new output directory
# ------------------------------------------------------------

dir.create(
  "08_results_primary",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 13. Save complete DESeq2 results
# ------------------------------------------------------------

write.table(
  res_df,
  file = "08_results_primary/deseq2_pupa_vs_larva_all_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Apply PRIMARY DEG criterion
# ------------------------------------------------------------

# padj < 0.01
# AND
# absolute log2FoldChange >= 1

deg_primary <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.01 &
    abs(res_df$log2FoldChange) >= 1,
]


# ------------------------------------------------------------
# 15. Separate pupa-up and larva-up DEGs
# ------------------------------------------------------------

pupa_up <- deg_primary[
  deg_primary$log2FoldChange >= 1,
]

larva_up <- deg_primary[
  deg_primary$log2FoldChange <= -1,
]


# ------------------------------------------------------------
# 16. Save DEG tables
# ------------------------------------------------------------

write.table(
  deg_primary,
  file = "08_results_primary/deseq2_pupa_vs_larva_DEGs_padj001_absLFC1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pupa_up,
  file = "08_results_primary/DEGs_pupa_up.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  larva_up,
  file = "08_results_primary/DEGs_larva_up.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Save normalized counts
# ------------------------------------------------------------

norm_counts <- counts(dds, normalized = TRUE)

norm_counts_df <- as.data.frame(norm_counts)

norm_counts_df$Geneid <- rownames(norm_counts_df)

norm_counts_df <- norm_counts_df[
  , c("Geneid", setdiff(colnames(norm_counts_df), "Geneid"))
]

write.table(
  norm_counts_df,
  file = "08_results_primary/normalized_counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Print DESeq2 summary
# ------------------------------------------------------------

cat("\n------------------------------\n")
cat("DESeq2 summary\n")
cat("------------------------------\n")

summary(res)


# ------------------------------------------------------------
# 19. Print primary DEG numbers
# ------------------------------------------------------------

cat("\n====================================\n")
cat("PRIMARY DEG CRITERION\n")
cat("padj < 0.01 and |log2FC| >= 1\n")
cat("====================================\n")

cat("Total primary DEGs :", nrow(deg_primary), "\n")
cat("Pupa-up DEGs       :", nrow(pupa_up), "\n")
cat("Larva-up DEGs      :", nrow(larva_up), "\n")
cat("Check total        :", nrow(pupa_up) + nrow(larva_up), "\n")