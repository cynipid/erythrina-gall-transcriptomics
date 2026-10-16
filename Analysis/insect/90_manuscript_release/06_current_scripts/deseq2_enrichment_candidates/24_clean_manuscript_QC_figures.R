# 24_clean_manuscript_QC_figures.R
# Purpose:
# Create manuscript-ready PCA and Pearson sample-correlation figures
# for the genome-guided Q. erythrinae larva-versus-pupa analysis.

library(DESeq2)
library(ggplot2)

project_dir <- "//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided"

count_file <- file.path(
  project_dir,
  "06_counts",
  "gene_counts_matrix.tsv"
)

metadata_file <- file.path(
  project_dir,
  "06_counts",
  "sample_metadata.tsv"
)

result_dir <- file.path(project_dir, "08_results_primary")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Read genome-guided gene counts and metadata
count_matrix <- read.delim(
  count_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

metadata <- read.delim(
  metadata_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Preserve the real sample identifiers used in the final analysis.
# Manuscript labels are separate so the source-library identities remain traceable.
sample_order <- c("L1", "L3", "L4", "P3", "P4", "P5")

sample_key <- data.frame(
  sample_id = sample_order,
  display_label = c("L1", "L2", "L3", "P1", "P2", "P3"),
  stringsAsFactors = FALSE
)

display_order <- sample_key$display_label

if (!setequal(colnames(count_matrix), sample_order)) {
  stop(
    "Count-matrix samples differ from the expected final samples: ",
    paste(colnames(count_matrix), collapse = ", ")
  )
}

if (!setequal(rownames(metadata), sample_order)) {
  stop(
    "Metadata samples differ from the expected final samples: ",
    paste(rownames(metadata), collapse = ", ")
  )
}

count_matrix <- count_matrix[, sample_order, drop = FALSE]
metadata <- metadata[sample_order, , drop = FALSE]

if (!identical(colnames(count_matrix), rownames(metadata))) {
  stop("Count-matrix columns and metadata rows are not aligned.")
}

metadata$stage <- tolower(trimws(as.character(metadata$stage)))
metadata$stage <- factor(
  metadata$stage,
  levels = c("larva", "pupa")
)

if (any(is.na(metadata$stage))) {
  stop("Metadata stage values must be larva or pupa.")
}

# Recreate the same filtered DESeq2 dataset used in the final analysis
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ stage
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Variance-stabilized expression values for sample-level QC
vsd <- vst(dds, blind = FALSE)

saveRDS(
  vsd,
  file.path(result_dir, "insect_gene_vsd.rds")
)

# ------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = "stage",
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_data$sample_id <- rownames(pca_data)

sample_key_index <- match(
  pca_data$sample_id,
  sample_key$sample_id
)

if (any(is.na(sample_key_index))) {
  stop("A PCA sample is missing from sample_key.")
}

pca_data$display_label <-
  sample_key$display_label[sample_key_index]

pca_data$stage <- factor(
  pca_data$stage,
  levels = c("larva", "pupa")
)

# Manual offsets keep labels readable for the tightly clustered samples
label_offsets <- data.frame(
  sample_id = sample_order,
  offset_x = c(3, -3, -3, -4, 4, 4),
  offset_y = c(-1, -2, 2, -2, 2, -2),
  stringsAsFactors = FALSE
)

offset_index <- match(
  pca_data$sample_id,
  label_offsets$sample_id
)

if (any(is.na(offset_index))) {
  stop("A PCA sample is missing from label_offsets.")
}

pca_data$label_PC1 <-
  pca_data$PC1 + label_offsets$offset_x[offset_index]

pca_data$label_PC2 <-
  pca_data$PC2 + label_offsets$offset_y[offset_index]

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = stage,
    shape = stage
  )
) +
  geom_point(size = 4) +
  geom_text(
    aes(
      x = label_PC1,
      y = label_PC2,
      label = display_label
    ),
    colour = "black",
    size = 4
  ) +
  scale_colour_manual(
    values = c(
      "larva" = "#009E73",
      "pupa" = "#CC79A7"
    ),
    labels = c(
      "larva" = "Larva",
      "pupa" = "Pupa"
    )
  ) +
  scale_shape_manual(
    values = c(
      "larva" = 16,
      "pupa" = 17
    ),
    labels = c(
      "larva" = "Larva",
      "pupa" = "Pupa"
    )
  ) +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    colour = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    plot.margin = margin(10, 25, 10, 10)
  )

write.table(
  pca_data,
  file.path(result_dir, "insect_gene_PCA_coordinates_clean.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_key,
  file.path(result_dir, "insect_sample_label_key.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  file.path(result_dir, "insect_gene_PCA_clean.png"),
  pca_plot,
  width = 8,
  height = 6,
  dpi = 400
)

ggsave(
  file.path(result_dir, "insect_gene_PCA_clean.pdf"),
  pca_plot,
  width = 8,
  height = 6
)

# ------------------------------------------------------------------
# Pearson sample-correlation heatmap
# ------------------------------------------------------------------

cor_mat <- cor(
  assay(vsd),
  method = "pearson"
)

cor_df <- as.data.frame(
  as.table(cor_mat),
  stringsAsFactors = FALSE
)

colnames(cor_df) <- c(
  "sample_row",
  "sample_column",
  "correlation"
)

cor_df$row_label <- sample_key$display_label[
  match(cor_df$sample_row, sample_key$sample_id)
]

cor_df$column_label <- sample_key$display_label[
  match(cor_df$sample_column, sample_key$sample_id)
]

if (any(is.na(cor_df$row_label)) || any(is.na(cor_df$column_label))) {
  stop("A correlation-matrix sample is missing from sample_key.")
}

cor_df$row_label <- factor(
  cor_df$row_label,
  levels = rev(display_order)
)

cor_df$column_label <- factor(
  cor_df$column_label,
  levels = display_order
)

cor_plot <- ggplot(
  cor_df,
  aes(
    x = column_label,
    y = row_label,
    fill = correlation
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = sprintf("%.2f", correlation)),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0.80,
    limits = c(0.60, 1.00),
    name = "Pearson r"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_fixed() +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

cor_mat_display <- cor_mat[sample_order, sample_order, drop = FALSE]
rownames(cor_mat_display) <- display_order
colnames(cor_mat_display) <- display_order

write.table(
  cor_mat_display,
  file.path(
    result_dir,
    "insect_gene_sample_correlation_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  cor_mat[sample_order, sample_order, drop = FALSE],
  file.path(
    result_dir,
    "insect_gene_sample_correlation_matrix_original_IDs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

ggsave(
  file.path(
    result_dir,
    "insect_gene_sample_correlation_heatmap.png"
  ),
  cor_plot,
  width = 6.5,
  height = 5.5,
  dpi = 400
)

ggsave(
  file.path(
    result_dir,
    "insect_gene_sample_correlation_heatmap.pdf"
  ),
  cor_plot,
  width = 6.5,
  height = 5.5
)

message("Clean insect PCA and correlation figures created.")

print(
  pca_data[, c(
    "sample_id",
    "display_label",
    "stage",
    "PC1",
    "PC2"
  )]
)

print(round(cor_mat_display, 3))
