# 06_add_GO_names_primary.R
# Add readable GO names and ontology classes
# to the significant primary DEG GO enrichment results.

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(GO.db)
library(AnnotationDbi)


# ------------------------------------------------------------
# 2. Create output folder
# ------------------------------------------------------------

dir.create(
  "10_functional_annotation/enrichment_results_primary_named",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Read significant GO enrichment tables
# ------------------------------------------------------------

GO_pupa <- read.delim(
  "10_functional_annotation/enrichment_results_primary/GO_enrichment_pupa_up_FDR005.tsv",
  header = TRUE,
  check.names = FALSE
)

GO_larva <- read.delim(
  "10_functional_annotation/enrichment_results_primary/GO_enrichment_larva_up_FDR005.tsv",
  header = TRUE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Collect all GO IDs
# ------------------------------------------------------------

all_GO_terms <- unique(
  c(
    GO_pupa$term,
    GO_larva$term
  )
)


# ------------------------------------------------------------
# 5. Retrieve GO names and ontology
# ------------------------------------------------------------

GO_names <- AnnotationDbi::select(
  GO.db,
  keys = all_GO_terms,
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
)

colnames(GO_names) <- c(
  "term",
  "term_name",
  "ontology"
)


# ------------------------------------------------------------
# 6. Add names to enrichment tables
# ------------------------------------------------------------

GO_pupa_named <- merge(
  GO_pupa,
  GO_names,
  by = "term",
  all.x = TRUE
)

GO_larva_named <- merge(
  GO_larva,
  GO_names,
  by = "term",
  all.x = TRUE
)


# Re-sort by significance

GO_pupa_named <- GO_pupa_named[
  order(
    GO_pupa_named$padj,
    GO_pupa_named$pvalue
  ),
]

GO_larva_named <- GO_larva_named[
  order(
    GO_larva_named$padj,
    GO_larva_named$pvalue
  ),
]


# ------------------------------------------------------------
# 7. Save named tables
# ------------------------------------------------------------

write.table(
  GO_pupa_named,
  "10_functional_annotation/enrichment_results_primary_named/GO_enrichment_pupa_up_FDR005_named.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  GO_larva_named,
  "10_functional_annotation/enrichment_results_primary_named/GO_enrichment_larva_up_FDR005_named.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
