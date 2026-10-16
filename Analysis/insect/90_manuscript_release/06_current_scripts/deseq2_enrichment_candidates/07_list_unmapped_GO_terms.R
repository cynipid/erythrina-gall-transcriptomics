setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")

pupa <- read.delim(
  "10_functional_annotation/enrichment_results_primary_named/GO_enrichment_pupa_up_FDR005_named.tsv",
  header = TRUE,
  check.names = FALSE
)

larva <- read.delim(
  "10_functional_annotation/enrichment_results_primary_named/GO_enrichment_larva_up_FDR005_named.tsv",
  header = TRUE,
  check.names = FALSE
)

pupa_unmapped <- pupa[is.na(pupa$term_name), ]
larva_unmapped <- larva[is.na(larva$term_name), ]

write.table(
  pupa_unmapped,
  "10_functional_annotation/enrichment_results_primary_named/GO_unmapped_pupa.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  larva_unmapped,
  "10_functional_annotation/enrichment_results_primary_named/GO_unmapped_larva.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

unique_unmapped <- unique(c(
  pupa_unmapped$term,
  larva_unmapped$term
))

write.table(
  data.frame(term = unique_unmapped),
  "10_functional_annotation/enrichment_results_primary_named/GO_unmapped_unique.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\nUnmapped GO terms\n")
cat("Pupa-up :", nrow(pupa_unmapped), "\n")
cat("Larva-up:", nrow(larva_unmapped), "\n")
cat("Unique  :", length(unique_unmapped), "\n")

head(unique_unmapped, 20)