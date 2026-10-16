# 09_volcano_plot_primary.R
# Volcano plot for official insect DESeq2 analysis
#
# Comparison: pupa vs larva
#
# Primary DEG definition:
# padj < 0.01 and |log2FC| >= 1
#
# Positive log2FC = Pupa-up
# Negative log2FC = Larva-up

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")

library(ggplot2)


# ------------------------------------------------------------
# 1. Read complete DESeq2 results
# ------------------------------------------------------------

res <- read.delim(
  "08_results_primary/deseq2_pupa_vs_larva_all_genes.tsv",
  header = TRUE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 2. Remove genes without adjusted P-values
# ------------------------------------------------------------

res_clean <- res[
  !is.na(res$padj) &
    !is.na(res$log2FoldChange),
]


# ------------------------------------------------------------
# 3. Classify genes
# ------------------------------------------------------------

res_clean$category <- "Not significant"

res_clean$category[
  res_clean$padj < 0.01 &
    res_clean$log2FoldChange >= 1
] <- "Pupa-up"

res_clean$category[
  res_clean$padj < 0.01 &
    res_clean$log2FoldChange <= -1
] <- "Larva-up"


# Fix legend order

res_clean$category <- factor(
  res_clean$category,
  levels = c(
    "Larva-up",
    "Pupa-up",
    "Not significant"
  )
)


# ------------------------------------------------------------
# 4. Calculate -log10 adjusted P-value
# ------------------------------------------------------------

res_clean$minus_log10_padj <- -log10(res_clean$padj)


# Cap extremely significant values at 50 for visualization
# This matches the style of the previous volcano figure.

res_clean$minus_log10_padj_plot <- pmin(
  res_clean$minus_log10_padj,
  50
)


# ------------------------------------------------------------
# 5. Create output folder
# ------------------------------------------------------------

dir.create(
  "09_figures_primary",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 6. Volcano plot
# ------------------------------------------------------------

volcano_plot <- ggplot(
  res_clean,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj_plot,
    color = category
  )
) +
  
  geom_point(
    alpha = 0.6,
    size = 1
  ) +
  
  scale_color_manual(
    values = c(
      "Larva-up" = "#009E73",
      "Pupa-up" = "#CC79A7",
      "Not significant" = "#BDBDBD"
    )
  ) +
  
  # log2FC thresholds: -1 and +1
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  
  # padj threshold: 0.01
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  
  scale_y_continuous(
    breaks = c(0, 10, 20, 30, 40, 50),
    labels = c("0", "10", "20", "30", "40", ">=50")
  ) +
  
  labs(
    x = expression(log[2]~"fold change (Pupa / Larva)"),
    y = expression(-log[10]~"(adjusted P value)"),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "right",
    
    axis.title = element_text(
      size = 14
    ),
    
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5,
      color = "black"
    ),
    
    legend.text = element_text(
      size = 12
    )
  )


# ------------------------------------------------------------
# 7. Save figure
# ------------------------------------------------------------

ggsave(
  "09_figures_primary/volcano_plot_pupa_vs_larva.png",
  plot = volcano_plot,
  width = 7,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  "09_figures_primary/volcano_plot_pupa_vs_larva.pdf",
  plot = volcano_plot,
  width = 7,
  height = 6,
  bg = "white"
)


# ------------------------------------------------------------
# 8. Short result summary
# ------------------------------------------------------------

table(res_clean$category)


# ------------------------------------------------------------
# 9. Show plot
# ------------------------------------------------------------

print(volcano_plot)