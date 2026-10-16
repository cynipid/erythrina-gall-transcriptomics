# 11_add_KEGG_names_primary.R
# Add readable KEGG pathway names to significant
# primary insect KEGG enrichment results.

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")

library(KEGGREST)


# ------------------------------------------------------------
# 1. Create output folder
# ------------------------------------------------------------

dir.create(
  "10_functional_annotation/enrichment_results_primary_named",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Read significant KEGG pathway tables
# ------------------------------------------------------------

KEGG_pupa <- read.delim(
  "10_functional_annotation/enrichment_results_primary/KEGG_pathway_enrichment_pupa_up_FDR005.tsv",
  header = TRUE,
  check.names = FALSE
)

KEGG_larva <- read.delim(
  "10_functional_annotation/enrichment_results_primary/KEGG_pathway_enrichment_larva_up_FDR005.tsv",
  header = TRUE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Retrieve KEGG pathway names
# ------------------------------------------------------------

pathway_names_vector <- KEGGREST::keggList("pathway")

pathway_names <- data.frame(
  term = sub(
    "^path:",
    "",
    names(pathway_names_vector)
  ),
  term_name = as.character(pathway_names_vector),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Add pathway names
# ------------------------------------------------------------

KEGG_pupa_named <- merge(
  KEGG_pupa,
  pathway_names,
  by = "term",
  all.x = TRUE
)

KEGG_larva_named <- merge(
  KEGG_larva,
  pathway_names,
  by = "term",
  all.x = TRUE
)


# Re-sort by significance

KEGG_pupa_named <- KEGG_pupa_named[
  order(
    KEGG_pupa_named$padj,
    KEGG_pupa_named$pvalue
  ),
]

KEGG_larva_named <- KEGG_larva_named[
  order(
    KEGG_larva_named$padj,
    KEGG_larva_named$pvalue
  ),
]


# ------------------------------------------------------------
# 5. Save named pathway tables
# ------------------------------------------------------------

write.table(
  KEGG_pupa_named,
  "10_functional_annotation/enrichment_results_primary_named/KEGG_pathway_enrichment_pupa_up_FDR005_named.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  KEGG_larva_named,
  "10_functional_annotation/enrichment_results_primary_named/KEGG_pathway_enrichment_larva_up_FDR005_named.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Short summary
# ------------------------------------------------------------

cat("\nKEGG naming summary\n")

cat(
  "Pupa-up pathways:",
  nrow(KEGG_pupa_named),
  "| unnamed:",
  sum(is.na(KEGG_pupa_named$term_name)),
  "\n"
)

cat(
  "Larva-up pathways:",
  nrow(KEGG_larva_named),
  "| unnamed:",
  sum(is.na(KEGG_larva_named$term_name)),
  "\n"
)

head(
  KEGG_pupa_named[
    , c("term", "term_name", "overlap", "padj")
  ],
  10
)

head(
  KEGG_larva_named[
    , c("term", "term_name", "overlap", "padj")
  ],
  10
)