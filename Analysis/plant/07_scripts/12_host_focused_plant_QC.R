#!/usr/bin/env Rscript
# Sample-level QC for the host-focused plant DESeq2 sensitivity branch.

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})

projects_root <- Sys.getenv("QE_PROJECTS_ROOT", unset = "")
if (!nzchar(projects_root)) {
  projects_root <- if (.Platform$OS.type == "windows") {
    "//wsl.localhost/Ubuntu/home/dodi/projects"
  } else file.path(Sys.getenv("HOME"), "projects")
}
root <- file.path(projects_root, "qe_cross_species_discovery",
                  "05_host_focused_DESeq2_sensitivity")
dds_file <- file.path(root, "results", "plant_gene_dds_host_focused.rds")
out_dir <- file.path(root, "qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stop_if <- function(x, message) if (x) stop(message, call. = FALSE)
stop_if(!file.exists(dds_file), paste("Missing host-focused DDS:", dds_file))

outs <- c(
  coordinates = file.path(out_dir, "plant_host_focused_PCA_coordinates.tsv"),
  correlation = file.path(out_dir, "plant_host_focused_sample_correlation.tsv"),
  pca_pdf = file.path(out_dir, "plant_host_focused_PCA.pdf"),
  correlation_pdf = file.path(out_dir, "plant_host_focused_sample_correlation_heatmap.pdf"),
  audit = file.path(out_dir, "plant_host_focused_QC_audit.txt")
)
dds <- readRDS(dds_file)
stop_if(!inherits(dds, "DESeqDataSet"), "Serialized object is not a DESeqDataSet")
stop_if(ncol(dds) != 6L, "Expected six plant samples")
expected <- c("LTS20_VA01_Sample_S7", "LTS20_VA01_Sample_S8",
              "LTS20_VA01_Sample_S9", "LTS20_VA01_Sample_S10",
              "LTS20_VA01_Sample_S11", "LTS20_VA01_Sample_S12")
stop_if(!identical(colnames(dds), expected), "Unexpected host-focused sample order")
stop_if(!identical(as.character(colData(dds)$condition), c("SG","SG","SG","S","S","S")),
        "Unexpected frozen condition order in host-focused DDS")

# Blind transformation is used solely for unsupervised PCA/correlation QC.
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca, "percentVar"), 2)
display <- c("GP1","GP2","GP3","UP1","UP2","UP3")
group <- factor(c("GP","GP","GP","UP","UP","UP"), levels = c("GP","UP"))
pca$sample_id <- rownames(pca)
pca$display_label <- display
pca$group <- group
write.table(pca, outs["coordinates"], sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

corr <- cor(assay(vsd), method = "pearson")
colnames(corr) <- rownames(corr) <- display
write.table(corr, outs["correlation"], sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)

pdf(outs["pca_pdf"], width = 7, height = 5.5, useDingbats = FALSE)
x_limits <- extendrange(pca$PC1, f = 0.12)
y_limits <- extendrange(pca$PC2, f = 0.16)
plot(pca$PC1, pca$PC2, col = c("#0072B2", "#D55E00")[group],
     pch = c(17,16)[group], cex = 1.5,
     xlim = x_limits, ylim = y_limits,
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"),
     main = "Host-focused plant transcriptome PCA")
text(pca$PC1, pca$PC2, labels = pca$display_label, pos = 3, cex = 0.9)
legend("topright", legend = c("Galled petiole", "Ungalled petiole"),
       col = c("#0072B2", "#D55E00"), pch = c(17,16), bty = "n")
dev.off()

pdf(outs["correlation_pdf"], width = 6.4, height = 5.8, useDingbats = FALSE)
par(mar = c(4,4,2.5,5))
image(1:6, 1:6, t(corr[6:1,]), col = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
      axes = FALSE, xlab = "", ylab = "", main = "Host-focused sample Pearson correlation")
axis(1, at = 1:6, labels = display)
axis(2, at = 1:6, labels = rev(display), las = 1)
for (i in 1:6) for (j in 1:6)
  text(i, 7-j, sprintf("%.2f", corr[j,i]), cex = 0.75)
box()
dev.off()

report <- c(
  "===== HOST-FOCUSED PLANT SAMPLE QC =====",
  paste("Host-focused genes:", nrow(dds)),
  "Transformation: DESeq2 vst(blind = TRUE) for unsupervised QC only.",
  "PCA: DESeq2 plotPCA with condition retained as display metadata.",
  "Correlation: Pearson correlation across VST-transformed gene profiles.",
  paste("PC1 variance:", percent_var[1], "%"),
  paste("PC2 variance:", percent_var[2], "%"),
  "Errors detected: 0",
  "HOST-FOCUSED PLANT SAMPLE QC: PASS"
)
writeLines(report, outs["audit"])
cat(paste(report, collapse = "\n"), "\n")
for (x in outs) cat("Output:", x, "\n")
