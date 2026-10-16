#!/usr/bin/env Rscript
# Narrow reviewer-requested QC sensitivity check for the frozen insect analysis.
# This script does not run DESeq(), change DE classifications, or overwrite
# the primary blind=FALSE figures and tables.
#
# In RStudio, open this file and click Source. The project is detected from the
# original Windows-to-WSL location or from the Linux WSL path. Optionally set
# EGW_PROJECT_DIR in the RStudio console before sourcing to override it.
# Command-line use is also supported:
#   Rscript insect_blind_TRUE_QC_sensitivity_RStudio.R /home/dodi/projects/qe_genome_guided

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) {
  project_dir <- args[[1L]]
} else {
  candidates <- unique(c(
    Sys.getenv("EGW_PROJECT_DIR", unset = ""),
    "//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided",
    "/home/dodi/projects/qe_genome_guided",
    getwd(), dirname(getwd())
  ))
  candidates <- candidates[nzchar(candidates)]
  available <- candidates[file.exists(file.path(
    candidates, "06_counts", "gene_counts_matrix.tsv"))]
  if (length(available) == 0L) {
    stop("Could not locate the frozen project inputs. In RStudio, set ",
         "EGW_PROJECT_DIR to the Q. erythrinae project root, then Source again.")
  }
  project_dir <- available[[1L]]
}
project_dir <- normalizePath(project_dir, mustWork = TRUE)

count_file <- file.path(project_dir, "06_counts", "gene_counts_matrix.tsv")
metadata_file <- file.path(project_dir, "06_counts", "sample_metadata.tsv")
primary_vsd_file <- file.path(project_dir, "08_results_primary", "insect_gene_vsd.rds")
out_dir <- file.path(project_dir, "review_sensitivity_20261011", "insect_blind_TRUE_QC")

for (path in c(count_file, metadata_file, primary_vsd_file)) {
  if (!file.exists(path)) stop("Missing frozen input: ", path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_order <- c("L1", "L3", "L4", "P3", "P4", "P5")
display_labels <- setNames(c("L1", "L2", "L3", "P1", "P2", "P3"), sample_order)
sample_stage <- setNames(c(rep("larva", 3L), rep("pupa", 3L)), sample_order)

count_data <- read.delim(count_file, header = TRUE, row.names = 1,
                         check.names = FALSE)
metadata <- read.delim(metadata_file, header = TRUE, row.names = 1,
                       check.names = FALSE)
if (!setequal(colnames(count_data), sample_order)) {
  stop("Count-matrix columns do not match the six frozen insect sample IDs.")
}
if (!setequal(rownames(metadata), sample_order) || !"stage" %in% names(metadata)) {
  stop("Sample metadata do not match the six frozen samples/stage design.")
}
count_data <- as.matrix(count_data[, sample_order, drop = FALSE])
metadata <- metadata[sample_order, , drop = FALSE]
if (anyNA(count_data) || any(count_data < 0) ||
    any(count_data != round(count_data))) {
  stop("The frozen count matrix must contain nonnegative integer counts.")
}
storage.mode(count_data) <- "integer"
metadata$stage <- factor(tolower(trimws(as.character(metadata$stage))),
                         levels = c("larva", "pupa"))
if (anyNA(metadata$stage) ||
    !identical(as.character(metadata$stage), unname(sample_stage))) {
  stop("Metadata stage labels do not match the frozen sample assignment.")
}

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata,
                              design = ~ stage)
dds <- dds[rowSums(counts(dds)) >= 10, ]
if (nrow(dds) < 1000L) stop("Fewer than 1000 genes passed the frozen filter.")

# The only new transformation in this review check is blind=TRUE.
vsd_true <- vst(dds, blind = TRUE)

summarize_qc <- function(vsd, mode) {
  mat <- assay(vsd)
  if (!identical(colnames(mat), sample_order)) {
    stop("Unexpected sample order in ", mode, " transformed matrix.")
  }
  pca <- plotPCA(vsd, intgroup = "stage", ntop = 500,
                 returnData = TRUE)
  pct <- 100 * attr(pca, "percentVar")
  pca$sample_id <- rownames(pca)
  pca$display_label <- unname(display_labels[pca$sample_id])
  if (anyNA(pca$display_label)) stop("Unrecognized PCA sample ID.")
  cm <- cor(mat, method = "pearson")
  pairs <- t(combn(sample_order, 2L))
  values <- cm[cbind(pairs[, 1L], pairs[, 2L])]
  within <- values[sample_stage[pairs[, 1L]] == sample_stage[pairs[, 2L]]]
  between <- values[sample_stage[pairs[, 1L]] != sample_stage[pairs[, 2L]]]
  lpc <- pca$PC1[pca$stage == "larva"]
  ppc <- pca$PC1[pca$stage == "pupa"]
  nonoverlap <- max(lpc) < min(ppc) || max(ppc) < min(lpc)
  metrics <- data.frame(
    transformation = mode,
    filtered_genes = nrow(mat),
    pca_ntop = 500L,
    PC1_variance_percent = pct[[1L]],
    PC2_variance_percent = pct[[2L]],
    larva_PC1_min = min(lpc), larva_PC1_max = max(lpc),
    pupa_PC1_min = min(ppc), pupa_PC1_max = max(ppc),
    PC1_ranges_nonoverlapping = nonoverlap,
    within_stage_r_min = min(within), within_stage_r_max = max(within),
    between_stage_r_min = min(between), between_stage_r_max = max(between),
    within_r_min_exceeds_between_r_max = min(within) > max(between)
  )
  list(pca = pca, cor = cm, metrics = metrics)
}

true_qc <- summarize_qc(vsd_true, "blind_TRUE")
vsd_false <- readRDS(primary_vsd_file)
if (!identical(rownames(assay(vsd_false)), rownames(assay(vsd_true)))) {
  stop("Primary blind=FALSE VST and blind=TRUE VST have different gene IDs.")
}
false_qc <- summarize_qc(vsd_false, "primary_blind_FALSE")
comparison <- rbind(false_qc$metrics, true_qc$metrics)
saveRDS(vsd_true, file.path(out_dir, "insect_vsd_blind_TRUE.rds"))

write.table(comparison,
            file.path(out_dir, "blind_TRUE_vs_primary_QC_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(true_qc$pca,
            file.path(out_dir, "insect_PCA_coordinates_blind_TRUE.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cm_display <- true_qc$cor[sample_order, sample_order, drop = FALSE]
rownames(cm_display) <- unname(display_labels[sample_order])
colnames(cm_display) <- unname(display_labels[sample_order])
write.table(cm_display,
            file.path(out_dir, "insect_correlation_matrix_blind_TRUE.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

p <- ggplot(true_qc$pca,
            aes(x = PC1, y = PC2, colour = stage, shape = stage)) +
  geom_point(size = 4) +
  geom_text(aes(label = display_label), hjust = -0.2, vjust = -0.5,
            colour = "black", size = 3.6, show.legend = FALSE) +
  scale_colour_manual(values = c(larva = "#009E73", pupa = "#CC79A7")) +
  scale_shape_manual(values = c(larva = 16, pupa = 17)) +
  labs(x = sprintf("PC1: %.1f%% variance", true_qc$metrics$PC1_variance_percent),
       y = sprintf("PC2: %.1f%% variance", true_qc$metrics$PC2_variance_percent),
       colour = NULL, shape = NULL) +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "insect_PCA_blind_TRUE.png"), p,
       width = 8, height = 6, dpi = 350)
ggsave(file.path(out_dir, "insect_PCA_blind_TRUE.pdf"), p,
       width = 8, height = 6)

cor_df <- as.data.frame(as.table(cm_display), stringsAsFactors = FALSE)
names(cor_df) <- c("row", "column", "r")
cor_df$row <- factor(cor_df$row, levels = rev(unname(display_labels)))
cor_df$column <- factor(cor_df$column, levels = unname(display_labels))
heat <- ggplot(cor_df, aes(x = column, y = row, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.7) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = mean(range(cor_df$r)),
                       limits = range(cor_df$r), name = "Pearson r") +
  coord_fixed() + labs(x = NULL, y = NULL) + theme_classic(base_size = 12) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())
ggsave(file.path(out_dir, "insect_correlation_blind_TRUE.png"), heat,
       width = 6.5, height = 5.5, dpi = 350)
ggsave(file.path(out_dir, "insect_correlation_blind_TRUE.pdf"), heat,
       width = 6.5, height = 5.5)

files_used <- c(count_file, metadata_file, primary_vsd_file)
write.table(data.frame(file = files_used,
                       md5 = unname(tools::md5sum(files_used))),
            file.path(out_dir, "input_md5.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "R_session_info.txt"))
print(comparison, row.names = FALSE)
message("Sensitivity check complete; outputs: ", out_dir)
