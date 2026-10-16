# 08_MA_plot_primary.R
# MA plot for official insect DESeq2 analysis
#
# DEG definition:
# padj < 0.01 and |log2FC| >= 1
#
# Color scheme:
# Larva-up        = green
# Pupa-up         = pink/purple
# Not significant = grey

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
# 2. Classify genes
# ------------------------------------------------------------

res$status <- "Not significant"

res$status[
  !is.na(res$padj) &
    res$padj < 0.01 &
    res$log2FoldChange >= 1
] <- "Pupa-up"

res$status[
  !is.na(res$padj) &
    res$padj < 0.01 &
    res$log2FoldChange <= -1
] <- "Larva-up"


# Fix order of categories in legend

res$status <- factor(
  res$status,
  levels = c(
    "Larva-up",
    "Pupa-up",
    "Not significant"
  )
)


# ------------------------------------------------------------
# 3. Prepare plotting data
# ------------------------------------------------------------

# Remove genes with baseMean = 0
# and genes without log2 fold-change estimates

plot_data <- res[
  !is.na(res$baseMean) &
    res$baseMean > 0 &
    !is.na(res$log2FoldChange),
]


# ------------------------------------------------------------
# 4. Create output folder
# ------------------------------------------------------------

dir.create(
  "09_figures_primary",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 5. Create MA plot
# ------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = baseMean,
    y = log2FoldChange,
    color = status
  )
) +
  
  geom_point(
    alpha = 0.5,
    size = 1
  ) +
  
  scale_color_manual(
    values = c(
      "Larva-up" = "#009E73",
      "Pupa-up" = "#CC79A7",
      "Not significant" = "#BDBDBD"
    )
  ) +
  
  scale_x_log10() +
  
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  
  labs(
    x = "Mean normalized expression",
    y = expression(log[2]~"fold change (Pupa / Larva)"),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.title = element_blank(),
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
# 6. Save figure
# ------------------------------------------------------------

ggsave(
  "09_figures_primary/MA_plot_pupa_vs_larva.png",
  plot = p,
  width = 7,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  "09_figures_primary/MA_plot_pupa_vs_larva.pdf",
  plot = p,
  width = 7,
  height = 6,
  bg = "white"
)


# ------------------------------------------------------------
# 7. Show plot in RStudio
# ------------------------------------------------------------

print(p)