# 03_check_primary_DEG_results.R
# Validate the official insect DEG table

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")

# Read primary DEG table
deg <- read.delim(
  "08_results_primary/deseq2_pupa_vs_larva_DEGs_padj001_absLFC1.tsv",
  header = TRUE,
  check.names = FALSE
)

cat("Total DEGs:", nrow(deg), "\n")

# Check adjusted P-value criterion
bad_padj <- deg[
  is.na(deg$padj) | deg$padj >= 0.01,
]

cat("Genes failing padj < 0.01:", nrow(bad_padj), "\n")

# Check fold-change criterion
bad_lfc <- deg[
  abs(deg$log2FoldChange) < 1,
]

cat("Genes failing |log2FC| >= 1:", nrow(bad_lfc), "\n")

# Count directions
pupa_up <- sum(deg$log2FoldChange >= 1)
larva_up <- sum(deg$log2FoldChange <= -1)

cat("Pupa-up:", pupa_up, "\n")
cat("Larva-up:", larva_up, "\n")
cat("Direction total:", pupa_up + larva_up, "\n")

# Final validation
if (
  nrow(bad_padj) == 0 &&
  nrow(bad_lfc) == 0 &&
  (pupa_up + larva_up) == nrow(deg)
) {
  cat("\nVALIDATION PASSED\n")
} else {
  cat("\nVALIDATION FAILED - inspect results\n")
}