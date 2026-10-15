#!/usr/bin/env Rscript

# Host-focused E. variegata GP-versus-UP MA and volcano plots.
# This is a visualization step only: it neither refits DESeq2 nor changes
# the frozen host-focused DEG calls.

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- "//wsl.localhost/Ubuntu/home/dodi/projects/qe_cross_species_discovery"
analysis_dir <- file.path(root, "05_host_focused_DESeq2_sensitivity")
input_file <- file.path(
  analysis_dir, "results",
  "plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv"
)
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Missing host-focused DESeq2 result table: ", input_file)

res <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("gene_id", "baseMean", "log2FoldChange", "padj")
missing <- setdiff(required, names(res))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
if (nrow(res) != 25779L) stop("Expected 25,779 host-focused result rows; found ", nrow(res))

res$status <- "Not differentially expressed"
res$status[is.na(res$padj)] <- "Not tested (NA padj)"
res$status[!is.na(res$padj) & res$padj < 0.01 & res$log2FoldChange >= 1] <- "GP-up"
res$status[!is.na(res$padj) & res$padj < 0.01 & res$log2FoldChange <= -1] <- "UP-up"
res$status <- factor(
  res$status,
  levels = c("GP-up", "UP-up", "Not differentially expressed", "Not tested (NA padj)")
)

expected_counts <- c("GP-up" = 3179L, "UP-up" = 1396L)
observed_counts <- table(res$status)
if (!all(observed_counts[names(expected_counts)] == expected_counts)) {
  stop("Host-focused DEG status counts do not match the frozen results.")
}

# Manuscript-wide condition colours: GP blue, UP orange.
status_colours <- c(
  "GP-up" = "#0072B2",
  "UP-up" = "#D55E00",
  "Not differentially expressed" = "#BDBDBD",
  "Not tested (NA padj)" = "#7A7A7A"
)

plot_data <- res[is.finite(res$baseMean) & is.finite(res$log2FoldChange), ]
plot_data$minus_log10_padj <- -log10(pmax(plot_data$padj, .Machine$double.xmin))
plot_data$minus_log10_padj[is.na(plot_data$padj)] <- NA_real_

common_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(0.75, "cm"),
    axis.title = element_text(face = "bold")
  )

legend_labels <- c(
  "GP-up" = "GP-up (n = 3,179)",
  "UP-up" = "UP-up (n = 1,396)",
  "Not differentially expressed" = "Not differentially expressed",
  "Not tested (NA padj)" = "Not tested"
)

ma_plot <- ggplot(plot_data, aes(x = log10(baseMean + 1), y = log2FoldChange, colour = status)) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"),
             colour = c("#777777", "#3D3D3D", "#777777"), linewidth = c(0.35, 0.4, 0.35)) +
  geom_point(data = subset(plot_data, status %in% c("Not differentially expressed", "Not tested (NA padj)")),
             size = 0.8, alpha = 0.34, na.rm = TRUE) +
  geom_point(data = subset(plot_data, status %in% c("GP-up", "UP-up")),
             size = 0.85, alpha = 0.56, na.rm = TRUE) +
  scale_colour_manual(values = status_colours, labels = legend_labels, drop = FALSE) +
  coord_cartesian(ylim = c(-12, 16)) +
  labs(
    x = expression(log[10](baseMean + 1)),
    y = expression(log[2]~fold~change~"(GP / UP)")
  ) +
  common_theme

volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = minus_log10_padj, colour = status)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "#777777", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "#777777", linewidth = 0.35) +
  geom_point(data = subset(plot_data, status %in% c("Not differentially expressed", "Not tested (NA padj)")),
             size = 0.8, alpha = 0.34, na.rm = TRUE) +
  geom_point(data = subset(plot_data, status %in% c("GP-up", "UP-up")),
             size = 0.85, alpha = 0.56, na.rm = TRUE) +
  scale_colour_manual(values = status_colours, labels = legend_labels, drop = FALSE) +
  labs(
    x = expression(log[2]~fold~change~"(GP / UP)"),
    y = expression(-log[10]~adjusted~P)
  ) +
  common_theme

save_plot <- function(plot, stem) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = 7.2, height = 5.4, units = "in")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = 7.2, height = 5.4,
         units = "in", dpi = 400)
}
save_plot(ma_plot, "plant_host_focused_MA_plot")
save_plot(volcano_plot, "plant_host_focused_volcano_plot")

status_output <- file.path(figure_dir, "plant_host_focused_DESeq2_plot_status_counts.tsv")
write.table(
  data.frame(status = names(observed_counts), genes = as.integer(observed_counts)),
  status_output, sep = "\t", quote = FALSE, row.names = FALSE
)

audit <- c(
  "===== HOST-FOCUSED MA AND VOLCANO PLOT BUILD =====",
  paste("Input:", input_file),
  paste("Result rows:", nrow(res)),
  paste("GP-up:", observed_counts[["GP-up"]]),
  paste("UP-up:", observed_counts[["UP-up"]]),
  "Colours: GP-up #0072B2; UP-up #D55E00; non-DE #BDBDBD; NA padj #7A7A7A",
  "No DESeq2 fitting or DEG classification was performed by this script.",
  "HOST-FOCUSED MA AND VOLCANO PLOT BUILD: PASS"
)
audit_output <- file.path(log_dir, "16_host_focused_MA_volcano_plot_audit.txt")
writeLines(audit, audit_output)
cat(paste(audit, collapse = "\n"), "\n")
cat("Output:", file.path(figure_dir, "plant_host_focused_MA_plot.pdf"), "\n")
cat("Output:", file.path(figure_dir, "plant_host_focused_volcano_plot.pdf"), "\n")
