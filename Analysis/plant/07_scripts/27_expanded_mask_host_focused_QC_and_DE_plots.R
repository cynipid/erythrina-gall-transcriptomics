#!/usr/bin/env Rscript

# Rebuild sample QC, MA and volcano figures for the expanded-mask host-focused
# plant analysis. This script reads the validated DESeq2 object/results and
# performs no model fitting, filtering, shrinkage or DEG reclassification.

options(stringsAsFactors = FALSE)

required_packages <- c("DESeq2", "SummarizedExperiment", "ggplot2")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
})

projects_root <- Sys.getenv("QE_PROJECTS_ROOT", unset = "")
if (!nzchar(projects_root)) {
  projects_root <- if (.Platform$OS.type == "windows") {
    "//wsl.localhost/Ubuntu/home/dodi/projects"
  } else {
    file.path(Sys.getenv("HOME"), "projects")
  }
}

discovery_root <- file.path(projects_root, "qe_cross_species_discovery")
branch_root <- file.path(discovery_root, "05_host_focused_DESeq2_sensitivity")
results_dir <- file.path(branch_root, "results")
figure_dir <- file.path(branch_root, "figures")
qc_dir <- file.path(branch_root, "qc")
logs_dir <- file.path(branch_root, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

dds_file <- file.path(
  results_dir, "plant_gene_dds_host_focused_expanded_mask.rds"
)
result_file <- file.path(
  results_dir,
  "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_all.tsv"
)

outputs <- list(
  pca_coordinates = file.path(
    qc_dir, "plant_host_focused_expanded_mask_PCA_coordinates.tsv"
  ),
  correlation = file.path(
    qc_dir, "plant_host_focused_expanded_mask_sample_correlation.tsv"
  ),
  pca_pdf = file.path(
    qc_dir, "plant_host_focused_expanded_mask_PCA.pdf"
  ),
  pca_png = file.path(
    qc_dir, "plant_host_focused_expanded_mask_PCA.png"
  ),
  correlation_pdf = file.path(
    qc_dir,
    "plant_host_focused_expanded_mask_sample_correlation_heatmap.pdf"
  ),
  correlation_png = file.path(
    qc_dir,
    "plant_host_focused_expanded_mask_sample_correlation_heatmap.png"
  ),
  ma_pdf = file.path(
    figure_dir, "plant_host_focused_expanded_mask_MA_plot.pdf"
  ),
  ma_png = file.path(
    figure_dir, "plant_host_focused_expanded_mask_MA_plot.png"
  ),
  volcano_pdf = file.path(
    figure_dir, "plant_host_focused_expanded_mask_volcano_plot.pdf"
  ),
  volcano_png = file.path(
    figure_dir, "plant_host_focused_expanded_mask_volcano_plot.png"
  ),
  status_counts = file.path(
    figure_dir,
    "plant_host_focused_expanded_mask_DESeq2_plot_status_counts.tsv"
  ),
  input_identity = file.path(
    logs_dir, "27_expanded_mask_host_focused_plot_input_identity.tsv"
  ),
  audit = file.path(
    logs_dir, "27_expanded_mask_host_focused_QC_and_DE_plot_audit.txt"
  )
)

audit_title <- "EXPANDED-MASK HOST-FOCUSED PLANT QC AND DE PLOT BUILD"
error_exit <- function(message) {
  cat(
    "===== ", audit_title, " =====\n",
    "ERROR: ", message, "\n",
    audit_title, ": CHECK REQUIRED\n",
    sep = "", file = stderr()
  )
  quit(status = 1L)
}

for (path in c(dds_file, result_file)) {
  if (!file.exists(path)) error_exit(paste("Missing required input:", path))
}
input_md5_before <- unname(tools::md5sum(c(dds_file, result_file)))

dds <- readRDS(dds_file)
if (!inherits(dds, "DESeqDataSet")) {
  error_exit("Serialized object is not a DESeqDataSet")
}
if (nrow(dds) != 24931L || ncol(dds) != 6L) {
  error_exit(paste(
    "Expected a 24931-gene by 6-sample DESeqDataSet; observed",
    nrow(dds), "by", ncol(dds)
  ))
}

expected_samples <- c(
  "LTS20_VA01_Sample_S7",
  "LTS20_VA01_Sample_S8",
  "LTS20_VA01_Sample_S9",
  "LTS20_VA01_Sample_S10",
  "LTS20_VA01_Sample_S11",
  "LTS20_VA01_Sample_S12"
)
display_labels <- c("GP1", "GP2", "GP3", "UP1", "UP2", "UP3")
display_groups <- factor(
  c("GP", "GP", "GP", "UP", "UP", "UP"),
  levels = c("GP", "UP")
)
if (!identical(colnames(dds), expected_samples)) {
  error_exit("DESeqDataSet sample order differs from the frozen six-sample order")
}
if (!identical(
      as.character(colData(dds)$condition),
      c("SG", "SG", "SG", "S", "S", "S")
    )) {
  error_exit("DESeqDataSet condition order differs from the frozen design")
}

res <- read.delim(
  result_file, sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("NA")
)
required_result_columns <- c(
  "gene_id", "baseMean", "log2FoldChange", "padj",
  "DESeq2_status_host_focused_expanded_mask"
)
missing_columns <- setdiff(required_result_columns, names(res))
if (length(missing_columns) > 0L) {
  error_exit(paste(
    "DESeq2 result table lacks columns:",
    paste(missing_columns, collapse = ", ")
  ))
}
if (nrow(res) != 24931L || anyDuplicated(res$gene_id) > 0L) {
  error_exit("Expanded-mask result table must contain 24931 unique genes")
}
if (!identical(res$gene_id, rownames(dds))) {
  error_exit("DESeq2 result and object gene orders are not identical")
}

calculated_status <- ifelse(
  is.na(res$padj),
  "NA_padj_not_testable",
  ifelse(
    res$padj < 0.01 & res$log2FoldChange >= 1,
    "GP_up",
    ifelse(
      res$padj < 0.01 & res$log2FoldChange <= -1,
      "UP_up",
      "non_DE"
    )
  )
)
if (!identical(calculated_status,
               res$DESeq2_status_host_focused_expanded_mask)) {
  error_exit("Stored DEG statuses disagree with the frozen reporting rule")
}
status_counts <- table(calculated_status)
expected_counts <- c(
  GP_up = 2349L,
  UP_up = 1386L,
  non_DE = 21022L,
  NA_padj_not_testable = 174L
)
for (status in names(expected_counts)) {
  observed <- if (status %in% names(status_counts)) {
    unname(status_counts[status])
  } else {
    0L
  }
  if (observed != expected_counts[status]) {
    error_exit(paste(
      "Unexpected", status, "count:", observed,
      "expected", expected_counts[status]
    ))
  }
}

# Blind VST is used only for unsupervised sample-level QC.
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(
  vsd, intgroup = "condition", ntop = min(500L, nrow(vsd)),
  returnData = TRUE
)
percent_var <- 100 * attr(pca, "percentVar")
if (nrow(pca) != 6L || length(percent_var) < 2L ||
    any(!is.finite(c(pca$PC1, pca$PC2, percent_var[1:2])))) {
  error_exit("PCA coordinates or explained-variance values are invalid")
}
if (!identical(rownames(pca), expected_samples)) {
  error_exit("PCA sample order differs from the frozen sample order")
}
pca$sample_id <- rownames(pca)
pca$display_label <- display_labels
pca$display_group <- display_groups
pca <- pca[, c(
  "sample_id", "display_label", "display_group", "condition", "PC1", "PC2"
)]
write.table(
  pca, outputs$pca_coordinates, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)

correlation <- cor(assay(vsd), method = "pearson")
if (!all(dim(correlation) == c(6L, 6L)) ||
    any(!is.finite(correlation)) ||
    any(abs(diag(correlation) - 1) > 1e-10)) {
  error_exit("VST sample-correlation matrix is invalid")
}
colnames(correlation) <- rownames(correlation) <- display_labels
write.table(
  correlation, outputs$correlation, sep = "\t", quote = FALSE,
  row.names = TRUE, col.names = NA
)

condition_colours <- c("GP" = "#0072B2", "UP" = "#D55E00")
condition_shapes <- c("GP" = 17, "UP" = 16)
status_levels <- c(
  "GP-up", "UP-up", "Not differentially expressed", "Not tested"
)
res$plot_status <- factor(
  ifelse(
    calculated_status == "GP_up", "GP-up",
    ifelse(
      calculated_status == "UP_up", "UP-up",
      ifelse(
        calculated_status == "NA_padj_not_testable",
        "Not tested", "Not differentially expressed"
      )
    )
  ),
  levels = status_levels
)
status_colours <- c(
  "GP-up" = "#0072B2",
  "UP-up" = "#D55E00",
  "Not differentially expressed" = "#BDBDBD",
  "Not tested" = "#7A7A7A"
)
legend_labels <- c(
  "GP-up" = "GP-up (n = 2,349)",
  "UP-up" = "UP-up (n = 1,386)",
  "Not differentially expressed" = "Not differentially expressed",
  "Not tested" = "Not tested"
)

common_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9.5),
    legend.key.width = grid::unit(0.7, "cm"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 16, 10, 10)
  )

pca_x_span <- max(diff(range(pca$PC1)), 1)
pca_y_span <- max(diff(range(pca$PC2)), 1)
pca$label_x <- pca$PC1
pca$label_y <- pca$PC2 + (0.045 * pca_y_span)

# UP replicates cluster tightly in PC space. Use deterministic offsets and
# leader lines to keep all three labels legible without moving any points.
up_idx <- which(pca$display_group == "UP")
if (length(up_idx) == 3L) {
  up_idx <- up_idx[order(pca$PC2[up_idx], decreasing = TRUE)]
  pca$label_x[up_idx] <- pca$PC1[up_idx] + c(0, -0.055, 0.055) * pca_x_span
  pca$label_y[up_idx] <- pca$PC2[up_idx] + c(0.080, 0.015, -0.080) * pca_y_span
}

pca_plot <- ggplot(
  pca, aes(x = PC1, y = PC2, colour = display_group, shape = display_group)
) +
  geom_hline(yintercept = 0, colour = "#E6E6E6", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "#E6E6E6", linewidth = 0.35) +
  geom_segment(
    aes(xend = label_x, yend = label_y),
    colour = "#777777", linewidth = 0.3, show.legend = FALSE
  ) +
  geom_point(size = 3.2) +
  geom_text(
    aes(x = label_x, y = label_y, label = display_label),
    show.legend = FALSE, size = 3.4, colour = "#222222"
  ) +
  scale_colour_manual(
    values = condition_colours,
    labels = c("Galled petiole (GP)", "Ungalled petiole (UP)"),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = condition_shapes,
    labels = c("Galled petiole (GP)", "Ungalled petiole (UP)"),
    drop = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.24))) +
  coord_cartesian(clip = "off") +
  labs(
    x = sprintf("PC1: %.2f%% variance", percent_var[1]),
    y = sprintf("PC2: %.2f%% variance", percent_var[2])
  ) +
  common_theme

correlation_long <- expand.grid(
  row_sample = display_labels,
  column_sample = display_labels,
  stringsAsFactors = FALSE
)
correlation_long$correlation <- mapply(
  function(row_name, column_name) correlation[row_name, column_name],
  correlation_long$row_sample,
  correlation_long$column_sample
)
correlation_long$row_sample <- factor(
  correlation_long$row_sample, levels = rev(display_labels)
)
correlation_long$column_sample <- factor(
  correlation_long$column_sample, levels = display_labels
)
correlation_limits <- range(correlation_long$correlation)
correlation_midpoint <- mean(correlation_limits)
correlation_plot <- ggplot(
  correlation_long,
  aes(x = column_sample, y = row_sample, fill = correlation)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.5) +
  scale_fill_gradient2(
    low = "#D55E00", mid = "#F7F7F7", high = "#0072B2",
    midpoint = correlation_midpoint,
    limits = correlation_limits,
    name = "Pearson r"
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

plot_data <- res[
  is.finite(res$baseMean) & is.finite(res$log2FoldChange), , drop = FALSE
]
plot_data$minus_log10_padj <- -log10(
  pmax(plot_data$padj, .Machine$double.xmin)
)
plot_data$minus_log10_padj[is.na(plot_data$padj)] <- NA_real_

ma_plot <- ggplot(
  plot_data,
  aes(x = log10(baseMean + 1), y = log2FoldChange, colour = plot_status)
) +
  geom_hline(
    yintercept = c(-1, 0, 1),
    linetype = c("dashed", "solid", "dashed"),
    colour = c("#777777", "#3D3D3D", "#777777"),
    linewidth = c(0.35, 0.4, 0.35)
  ) +
  geom_point(
    data = subset(
      plot_data,
      plot_status %in% c("Not differentially expressed", "Not tested")
    ),
    size = 0.8, alpha = 0.34, na.rm = TRUE
  ) +
  geom_point(
    data = subset(plot_data, plot_status %in% c("GP-up", "UP-up")),
    size = 0.85, alpha = 0.56, na.rm = TRUE
  ) +
  scale_colour_manual(
    values = status_colours, labels = legend_labels, drop = FALSE
  ) +
  labs(
    x = expression(log[10](baseMean + 1)),
    y = expression(log[2]~fold~change~"(GP / UP)")
  ) +
  common_theme

volcano_plot <- ggplot(
  plot_data,
  aes(x = log2FoldChange, y = minus_log10_padj, colour = plot_status)
) +
  geom_vline(
    xintercept = c(-1, 1), linetype = "dashed",
    colour = "#777777", linewidth = 0.35
  ) +
  geom_hline(
    yintercept = -log10(0.01), linetype = "dashed",
    colour = "#777777", linewidth = 0.35
  ) +
  geom_point(
    data = subset(
      plot_data,
      plot_status %in% c("Not differentially expressed", "Not tested")
    ),
    size = 0.8, alpha = 0.34, na.rm = TRUE
  ) +
  geom_point(
    data = subset(plot_data, plot_status %in% c("GP-up", "UP-up")),
    size = 0.85, alpha = 0.56, na.rm = TRUE
  ) +
  scale_colour_manual(
    values = status_colours, labels = legend_labels, drop = FALSE
  ) +
  labs(
    x = expression(log[2]~fold~change~"(GP / UP)"),
    y = expression(-log[10]~adjusted~P)
  ) +
  common_theme

save_plot <- function(plot_object, pdf_path, png_path, width, height) {
  ggsave(
    pdf_path, plot_object, width = width, height = height,
    units = "in", device = cairo_pdf
  )
  ggsave(
    png_path, plot_object, width = width, height = height,
    units = "in", dpi = 400, bg = "white"
  )
}
save_plot(pca_plot, outputs$pca_pdf, outputs$pca_png, 7.2, 5.4)
save_plot(
  correlation_plot, outputs$correlation_pdf, outputs$correlation_png,
  6.6, 5.8
)
save_plot(ma_plot, outputs$ma_pdf, outputs$ma_png, 7.2, 5.4)
save_plot(volcano_plot, outputs$volcano_pdf, outputs$volcano_png, 7.2, 5.4)

status_output <- data.frame(
  status = c("GP-up", "UP-up", "Not differentially expressed", "Not tested"),
  genes = c(
    expected_counts["GP_up"], expected_counts["UP_up"],
    expected_counts["non_DE"], expected_counts["NA_padj_not_testable"]
  ),
  stringsAsFactors = FALSE
)
write.table(
  status_output, outputs$status_counts, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)

input_identity <- data.frame(
  input = c("expanded_mask_dds", "expanded_mask_DESeq2_results"),
  path = c(dds_file, result_file),
  MD5 = input_md5_before,
  stringsAsFactors = FALSE
)
write.table(
  input_identity, outputs$input_identity, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)

input_md5_after <- unname(tools::md5sum(c(dds_file, result_file)))
if (!identical(input_md5_before, input_md5_after)) {
  error_exit("One or more plot inputs changed during figure generation")
}
missing_outputs <- unlist(outputs)[!file.exists(unlist(outputs))]
# The audit is written below, so it is the only permissible missing output here.
missing_before_audit <- setdiff(missing_outputs, outputs$audit)
if (length(missing_before_audit) > 0L) {
  error_exit(paste(
    "Expected output files were not generated:",
    paste(missing_before_audit, collapse = ", ")
  ))
}

audit_lines <- c(
  paste0("===== ", audit_title, " ====="),
  paste("Expanded-mask genes represented in DESeq2 object:", nrow(dds)),
  "Sample labels: S7->GP1, S8->GP2, S9->GP3, S10->UP1, S11->UP2, S12->UP3",
  "Legacy metadata translation: SG=GP and S=UP",
  "Transformation: DESeq2 vst(blind=TRUE), used only for unsupervised QC.",
  paste("PCA genes used:", min(500L, nrow(vsd))),
  sprintf("PC1 variance: %.2f%%", percent_var[1]),
  sprintf("PC2 variance: %.2f%%", percent_var[2]),
  "PCA labels: deterministic separated offsets with leader lines; sample coordinates unchanged.",
  "Correlation: Pearson correlation across all VST-transformed tested genes.",
  paste("GP-up:", expected_counts["GP_up"]),
  paste("UP-up:", expected_counts["UP_up"]),
  paste("Non-DE:", expected_counts["non_DE"]),
  paste("Not tested (NA padj):", expected_counts["NA_padj_not_testable"]),
  "Colours: GP #0072B2; UP #D55E00; non-DE #BDBDBD; not tested #7A7A7A",
  "No DESeq2 fitting, filtering, shrinkage or DEG reclassification was performed.",
  "Input MD5 verification after output generation: PASS",
  "Errors detected: 0",
  paste0(audit_title, ": PASS")
)
writeLines(audit_lines, outputs$audit)
cat(paste(audit_lines, collapse = "\n"), "\n")
for (path in unlist(outputs)) cat("Output:", path, "\n")
