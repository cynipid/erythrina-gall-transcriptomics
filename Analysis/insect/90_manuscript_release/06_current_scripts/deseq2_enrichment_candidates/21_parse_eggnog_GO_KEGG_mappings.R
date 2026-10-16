# 21_parse_eggnog_GO_KEGG_mappings.R
# Purpose:
# Parse eggNOG-mapper output and create gene-to-GO / KEGG mapping tables
# for enrichment analysis.

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")

dir.create("10_functional_annotation/parsed_mappings", showWarnings = FALSE)

eggnog_file <- "10_functional_annotation/eggnog_out/Qerythrinae.emapper.annotations"

# -----------------------------
# 1. Read eggNOG annotations
# -----------------------------

eggnog_lines <- readLines(eggnog_file, warn = FALSE)

# Find header line
header_line <- eggnog_lines[grepl("^#query", eggnog_lines)]

if (length(header_line) != 1) {
  stop("Could not find unique #query header line in eggNOG output.")
}

column_names <- strsplit(header_line, "\t")[[1]]
column_names[1] <- sub("^#", "", column_names[1])

# Data lines: remove all comment lines
data_lines <- eggnog_lines[!grepl("^#", eggnog_lines)]

eggnog <- read.delim(
  textConnection(data_lines),
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(eggnog) <- column_names

# Rename query column to Geneid for consistency with DESeq2
colnames(eggnog)[colnames(eggnog) == "query"] <- "Geneid"

# -----------------------------
# 2. Helper function for parsing comma-separated annotation terms
# -----------------------------

parse_terms <- function(data, term_column, output_term_name) {
  
  temp <- data[, c("Geneid", term_column)]
  colnames(temp) <- c("Geneid", "terms")
  
  # Remove missing terms
  temp <- temp[!is.na(temp$terms) & temp$terms != "" & temp$terms != "-", ]
  
  if (nrow(temp) == 0) {
    return(data.frame(Geneid = character(), term = character()))
  }
  
  split_terms <- strsplit(temp$terms, ",")
  
  parsed <- data.frame(
    Geneid = rep(temp$Geneid, lengths(split_terms)),
    term = unlist(split_terms),
    stringsAsFactors = FALSE
  )
  
  parsed$term <- trimws(parsed$term)
  parsed <- parsed[parsed$term != "" & parsed$term != "-", ]
  
  colnames(parsed) <- c("Geneid", output_term_name)
  
  # Remove duplicate gene-term pairs
  parsed <- unique(parsed)
  
  return(parsed)
}

# -----------------------------
# 3. Parse GO / KEGG mappings
# -----------------------------

gene_to_GO <- parse_terms(
  data = eggnog,
  term_column = "GOs",
  output_term_name = "GO"
)

gene_to_KEGG_ko <- parse_terms(
  data = eggnog,
  term_column = "KEGG_ko",
  output_term_name = "KEGG_ko"
)

gene_to_KEGG_pathway <- parse_terms(
  data = eggnog,
  term_column = "KEGG_Pathway",
  output_term_name = "KEGG_Pathway"
)

# -----------------------------
# 4. Save mapping files
# -----------------------------

write.table(
  eggnog,
  file = "10_functional_annotation/parsed_mappings/Qerythrinae_emapper_annotations_parsed.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_to_GO,
  file = "10_functional_annotation/parsed_mappings/gene_to_GO.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_to_KEGG_ko,
  file = "10_functional_annotation/parsed_mappings/gene_to_KEGG_ko.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_to_KEGG_pathway,
  file = "10_functional_annotation/parsed_mappings/gene_to_KEGG_pathway.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 5. Summary table
# -----------------------------

mapping_summary <- data.frame(
  item = c(
    "eggNOG annotation rows",
    "Unique annotated genes in eggNOG",
    "Genes with at least one GO term",
    "GO gene-term pairs",
    "Unique GO terms",
    "Genes with at least one KEGG KO",
    "KEGG KO gene-term pairs",
    "Unique KEGG KO terms",
    "Genes with at least one KEGG pathway",
    "KEGG pathway gene-term pairs",
    "Unique KEGG pathways"
  ),
  count = c(
    nrow(eggnog),
    length(unique(eggnog$Geneid)),
    length(unique(gene_to_GO$Geneid)),
    nrow(gene_to_GO),
    length(unique(gene_to_GO$GO)),
    length(unique(gene_to_KEGG_ko$Geneid)),
    nrow(gene_to_KEGG_ko),
    length(unique(gene_to_KEGG_ko$KEGG_ko)),
    length(unique(gene_to_KEGG_pathway$Geneid)),
    nrow(gene_to_KEGG_pathway),
    length(unique(gene_to_KEGG_pathway$KEGG_Pathway))
  )
)

write.table(
  mapping_summary,
  file = "10_functional_annotation/parsed_mappings/eggnog_mapping_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 6. Print checks
# -----------------------------

mapping_summary

head(gene_to_GO, 10)
head(gene_to_KEGG_ko, 10)
head(gene_to_KEGG_pathway, 10)

list.files("10_functional_annotation/parsed_mappings")