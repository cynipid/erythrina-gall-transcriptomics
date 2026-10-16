# ============================================================================
# MEROPS cross-species quantitative figure, curated v2
#
# Run in RStudio:
# source("//wsl.localhost/Ubuntu/home/dodi/projects/qe_merops_screen/scripts/13_merops_cross_species_quantitative_figure_v2.R")
#
# Outputs one four-panel manuscript figure plus separate panel files.
# The plots compare candidate programmes; they do not assert biochemical
# activity, host delivery, or direct plant--insect protein interactions.
# ============================================================================

options(stringsAsFactors = FALSE)

# 1. Paths -----------------------------------------------------------------

base_dir <- "//wsl.localhost/Ubuntu/home/dodi/projects/qe_merops_screen"
results_dir <- file.path(base_dir, "results")
figure_dir <- file.path(base_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

family_file <- file.path(
  results_dir,
  "MEROPS_insect_inhibitor_family_by_DESeq2_curated_v2.tsv"
)
larva_file <- file.path(
  results_dir,
  "insect_larva_up_inhibitor_annotation_domain_audit_v2.tsv"
)
plant_file <- file.path(
  results_dir,
  "plant_GPup_MEROPS_peptidase_priority_curated_v2.tsv"
)
fisher_file <- file.path(
  results_dir,
  "MEROPS_insect_secreted_inhibitor_direction_sensitivity_v2.tsv"
)

required_files <- c(family_file, larva_file, plant_file, fisher_file)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required input file(s):\n",
    paste(required_files[!file.exists(required_files)], collapse = "\n")
  )
}

# 2. Packages and helpers ---------------------------------------------------

required_packages <- c("ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the missing package(s), restart R, and source this script again: ",
    paste(missing_packages, collapse = ", "),
    "\nUse: install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "))"
  )
}

library(ggplot2)
library(patchwork)

read_tsv <- function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
}

require_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

as_numeric_safely <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

is_priority <- function(x) {
  grepl("^priority_[12]_", as.character(x))
}

plant_family_group <- function(x) {
  x <- as.character(x)
  ifelse(
    grepl("^A01", x), "A01: aspartic peptidase A1",
    ifelse(
      grepl("^C01", x), "C01: papain-like cysteine peptidase",
      ifelse(
        grepl("^S08", x), "S08: subtilisin-like serine peptidase",
        ifelse(
          grepl("^S10", x), "S10: serine carboxypeptidase",
          ifelse(
            grepl("^M10", x), "M10: metallopeptidase",
            ifelse(
              grepl("^C13", x), "C13: legumain/VPE",
              ifelse(
                grepl("^S28", x), "S28: serine carboxypeptidase",
                x
              )
            )
          )
        )
      )
    )
  )
}

theme_merops <- function(base_size = 10.5) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey94", colour = "grey70"),
      strip.text = element_text(face = "bold")
    )
}

# Colour choices are fixed to remain stable between reruns.
status_colours <- c(
  pupa_up = "#D55E00",
  larva_up = "#0072B2",
  non_DE = "#A8A8A8",
  NA_padj_not_testable = "#CC79A7",
  excluded_low_count = "#E6AB02"
)

larva_family_colours <- c(
  I01 = "#0072B2",
  I04 = "#D55E00",
  I19 = "#009E73"
)

# 3. Read and validate ------------------------------------------------------

family <- read_tsv(family_file)
larva <- read_tsv(larva_file)
plant_all <- read_tsv(plant_file)
fisher_all <- read_tsv(fisher_file)

require_columns(
  family,
  c("merops_family", "DESeq2_status", "candidate_genes"),
  "Insect family summary"
)
require_columns(
  larva,
  c("gene_id", "merops_family", "log2FoldChange", "padj"),
  "Larva-up annotation/domain table"
)
require_columns(
  plant_all,
  c(
    "gene_id", "merops_family", "manuscript_priority",
    "log2FoldChange", "padj"
  ),
  "Plant priority table"
)
require_columns(
  fisher_all,
  c(
    "candidate_set", "odds_ratio", "confidence_low", "confidence_high",
    "pvalue", "padj_BH_across_three_nested_tests"
  ),
  "Insect direction-sensitivity table"
)

family$candidate_genes <- as_numeric_safely(family$candidate_genes)
larva$log2FoldChange <- as_numeric_safely(larva$log2FoldChange)
larva$padj <- as_numeric_safely(larva$padj)
plant_all$log2FoldChange <- as_numeric_safely(plant_all$log2FoldChange)
plant_all$padj <- as_numeric_safely(plant_all$padj)
fisher_all$odds_ratio <- as_numeric_safely(fisher_all$odds_ratio)
fisher_all$confidence_low <- as_numeric_safely(fisher_all$confidence_low)
fisher_all$confidence_high <- as_numeric_safely(fisher_all$confidence_high)
fisher_all$pvalue <- as_numeric_safely(fisher_all$pvalue)
fisher_all$padj_BH_across_three_nested_tests <- as_numeric_safely(
  fisher_all$padj_BH_across_three_nested_tests
)

plant <- plant_all[is_priority(plant_all$manuscript_priority), , drop = FALSE]

errors <- character()
if (sum(family$candidate_genes, na.rm = TRUE) != 149L) {
  errors <- c(errors, "Insect family table does not total 149 candidates")
}
if (nrow(larva) != 9L || length(unique(larva$gene_id)) != 9L) {
  errors <- c(errors, "Larva table does not contain 9 unique candidates")
}
if (!setequal(unique(larva$merops_family), c("I01", "I04", "I19"))) {
  errors <- c(errors, "Larva table does not contain only I01, I04 and I19")
}
if (nrow(plant) != 20L || length(unique(plant$gene_id)) != 20L) {
  errors <- c(errors, "Plant priority subset does not contain 20 unique candidates")
}
if (any(!is.finite(c(larva$log2FoldChange, plant$log2FoldChange)))) {
  errors <- c(errors, "One or more candidate log2 fold changes is non-numeric")
}

if (length(errors)) {
  cat("\nPre-plot errors detected:", length(errors), "\n")
  for (error in errors) cat("ERROR:", error, "\n")
  stop("MEROPS quantitative figure input audit failed")
}

# 4A. All insect inhibitor families by DESeq2 status -----------------------

status_order <- c(
  "pupa_up",
  "larva_up",
  "non_DE",
  "NA_padj_not_testable",
  "excluded_low_count"
)
status_labels <- c(
  pupa_up = "Pupa-up",
  larva_up = "Larva-up",
  non_DE = "Not a primary DEG",
  NA_padj_not_testable = "NA adjusted P",
  excluded_low_count = "Excluded: low count"
)

family_totals <- aggregate(
  candidate_genes ~ merops_family,
  data = family,
  FUN = sum
)
family_levels <- family_totals$merops_family[
  order(family_totals$candidate_genes, family_totals$merops_family)
]

family$merops_family <- factor(family$merops_family, levels = family_levels)
family$DESeq2_status <- factor(family$DESeq2_status, levels = status_order)

panel_a <- ggplot(
  family,
  aes(x = merops_family, y = candidate_genes, fill = DESeq2_status)
) +
  geom_col(width = 0.78) +
  coord_flip() +
  scale_fill_manual(
    values = status_colours,
    breaks = status_order,
    labels = status_labels,
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = "Insect secreted inhibitor-family homologues",
    subtitle = "Curated A+B candidates (n = 149); stacked by DESeq2 status",
    x = "MEROPS inhibitor family",
    y = "Candidate genes",
    fill = "Expression status"
  ) +
  theme_merops() +
  theme(legend.position = "bottom")

# 4B. Focused insect larva-up candidates ----------------------------------

larva$larva_vs_pupa_log2FC <- -larva$log2FoldChange
larva$display_gene <- factor(
  larva$gene_id,
  levels = larva$gene_id[order(larva$larva_vs_pupa_log2FC)]
)

panel_b <- ggplot(
  larva,
  aes(
    x = larva_vs_pupa_log2FC,
    y = display_gene,
    colour = merops_family
  )
) +
  geom_segment(
    aes(x = 0, xend = larva_vs_pupa_log2FC, yend = display_gene),
    linewidth = 0.55,
    colour = "grey72"
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(label = sprintf("%.2f", larva_vs_pupa_log2FC)),
    hjust = -0.2,
    size = 3,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  scale_colour_manual(
    values = larva_family_colours,
    labels = c(
      I01 = "I01 Kazal-type",
      I04 = "I04 serpin",
      I19 = "I19 pacifastin"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.14))) +
  labs(
    title = "Focused larva-up inhibitor candidates",
    subtitle = "Nine candidates with the expected Pfam-family support",
    x = expression(log[2]~fold~change~"(larva / pupa)"),
    y = "Insect gene",
    colour = "MEROPS family"
  ) +
  theme_merops()

# 4C. Priority plant GP-up peptidase-related candidates --------------------

plant$family_group <- plant_family_group(plant$merops_family)
plant$display_gene <- factor(
  plant$gene_id,
  levels = plant$gene_id[order(plant$log2FoldChange)]
)

plant_family_levels <- sort(unique(plant$family_group))
plant_palette <- setNames(
  grDevices::hcl.colors(
    length(plant_family_levels),
    palette = "Dark 3",
    rev = FALSE
  ),
  plant_family_levels
)

panel_c <- ggplot(
  plant,
  aes(x = log2FoldChange, y = display_gene, colour = family_group)
) +
  geom_segment(
    aes(x = 0, xend = log2FoldChange, yend = display_gene),
    linewidth = 0.48,
    colour = "grey72"
  ) +
  geom_point(size = 2.7) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  scale_colour_manual(values = plant_palette) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  labs(
    title = "Galled-petiole-up plant peptidase candidates",
    subtitle = "Priority 1 + 2 classically secreted candidates (n = 20)",
    x = expression(log[2]~fold~change~"(galled / ungalled petiole)"),
    y = "Plant gene",
    colour = "MEROPS group"
  ) +
  theme_merops(9.5) +
  theme(
    axis.text.y = element_text(size = 7.8),
    legend.text = element_text(size = 8)
  )

# 4D. Forest plot of insect pupa-versus-larva enrichment ------------------

forest_sets <- c(
  "curated_A_plus_B_inhibitor_candidates",
  "Tier_A_named_reference_matches"
)
forest <- fisher_all[
  fisher_all$candidate_set %in% forest_sets,
  , drop = FALSE
]
forest$display_set <- ifelse(
  forest$candidate_set == "curated_A_plus_B_inhibitor_candidates",
  "All curated A+B candidates",
  "Tier A named-reference matches"
)
forest$display_set <- factor(
  forest$display_set,
  levels = c(
    "All curated A+B candidates",
    "Tier A named-reference matches"
  )
)
forest$fdr_label <- sprintf(
  "OR %.2f; FDR %.2g",
  forest$odds_ratio,
  forest$padj_BH_across_three_nested_tests
)

if (nrow(forest) != 2L ||
    any(!is.finite(c(
      forest$odds_ratio,
      forest$confidence_low,
      forest$confidence_high
    )))) {
  stop("Forest plot could not identify the two expected valid comparisons")
}

panel_d <- ggplot(
  forest,
  aes(x = display_set, y = odds_ratio)
) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_errorbar(
    aes(ymin = confidence_low, ymax = confidence_high),
    width = 0.15,
    linewidth = 0.7,
    colour = "#333333"
  ) +
  geom_point(size = 3.2, colour = "#7A3E9D") +
  geom_text(
    aes(y = 16.5, label = fdr_label),
    hjust = 1,
    nudge_x = 0.16,
    size = 3.1
  ) +
  coord_flip() +
  scale_y_log10(
    breaks = c(1, 2, 4, 8, 16),
    limits = c(0.75, 18)
  ) +
  labs(
    title = "Insect inhibitor-set stage bias",
    subtitle = "Fisher exact tests; values above 1 indicate pupa-up enrichment",
    x = NULL,
    y = "Odds ratio (log scale; 95% CI)"
  ) +
  theme_merops() +
  theme(legend.position = "none")

# 5. Save separate panels --------------------------------------------------

panel_files <- c(
  insect_family = file.path(figure_dir, "MEROPS_panel_A_insect_family_status_v2"),
  larva_candidates = file.path(figure_dir, "MEROPS_panel_B_larva_candidates_v2"),
  plant_candidates = file.path(figure_dir, "MEROPS_panel_C_plant_candidates_v2"),
  stage_bias = file.path(figure_dir, "MEROPS_panel_D_insect_stage_bias_v2")
)

save_panel <- function(plot_object, path_without_extension, width, height) {
  ggsave(
    paste0(path_without_extension, ".png"),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 450,
    bg = "white"
  )
  ggsave(
    paste0(path_without_extension, ".pdf"),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}

save_panel(panel_a, panel_files[["insect_family"]], 8.2, 7.2)
save_panel(panel_b, panel_files[["larva_candidates"]], 8.2, 5.8)
save_panel(panel_c, panel_files[["plant_candidates"]], 9.2, 8.3)
save_panel(panel_d, panel_files[["stage_bias"]], 8.2, 4.5)

# 6. Assemble and save the manuscript figure ------------------------------

top_row <- panel_a + panel_b + plot_layout(widths = c(1.05, 1))
bottom_row <- panel_c + panel_d + plot_layout(widths = c(1.18, 0.82))

combined <- (top_row / bottom_row) +
  plot_annotation(
    title = "Stage-associated insect inhibitor and plant peptidase programmes",
    subtitle = paste(
      "Homology- and domain-supported candidate sets;",
      "no direct plant--insect protein interaction is inferred"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, colour = "grey30"),
      plot.tag = element_text(face = "bold", size = 14)
    )
  ) &
  theme(plot.tag.position = c(0.01, 0.99))

combined_png <- file.path(
  figure_dir,
  "MEROPS_cross_species_quantitative_figure_curated_v2.png"
)
combined_pdf <- file.path(
  figure_dir,
  "MEROPS_cross_species_quantitative_figure_curated_v2.pdf"
)
audit_output <- file.path(
  figure_dir,
  "MEROPS_cross_species_quantitative_figure_curated_v2_audit.tsv"
)

ggsave(
  combined_png,
  plot = combined,
  width = 16,
  height = 13,
  units = "in",
  dpi = 450,
  bg = "white"
)
ggsave(
  combined_pdf,
  plot = combined,
  width = 16,
  height = 13,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

# 7. Final audit -----------------------------------------------------------

audit <- data.frame(
  metric = c(
    "insect_curated_candidates_in_family_table",
    "insect_unique_families",
    "insect_larva_up_Pfam_supported_candidates",
    "plant_GP_up_priority_candidates",
    "forest_plot_comparisons",
    "combined_PNG_exists",
    "combined_PDF_exists"
  ),
  observed = c(
    sum(family$candidate_genes),
    length(unique(family$merops_family)),
    nrow(larva),
    nrow(plant),
    nrow(forest),
    as.integer(file.exists(combined_png)),
    as.integer(file.exists(combined_pdf))
  ),
  expected = c(149L, 18L, 9L, 20L, 2L, 1L, 1L),
  stringsAsFactors = FALSE
)
audit$match <- audit$observed == audit$expected

write.table(
  audit,
  audit_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n===== MEROPS CROSS-SPECIES QUANTITATIVE FIGURE V2 =====\n")
print(audit, row.names = FALSE)
cat("\nCombined PNG:", combined_png, "\n")
cat("Combined PDF:", combined_pdf, "\n")
cat("Audit:", audit_output, "\n")

if (!all(audit$match)) {
  cat("\nMEROPS CROSS-SPECIES QUANTITATIVE FIGURE V2: CHECK REQUIRED\n")
  stop("One or more final figure audit checks failed")
}

cat("\nMEROPS CROSS-SPECIES QUANTITATIVE FIGURE V2: PASS\n")
