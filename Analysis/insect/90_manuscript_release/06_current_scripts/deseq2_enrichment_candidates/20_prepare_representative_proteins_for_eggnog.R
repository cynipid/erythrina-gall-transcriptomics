# 20_prepare_representative_proteins_for_eggnog.R
# Purpose:
# Prepare one representative protein per gene for eggNOG-mapper.
# Rule:
# For each gene/locus_tag, keep the longest protein isoform.

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")

dir.create("10_functional_annotation", showWarnings = FALSE)
dir.create("10_functional_annotation/eggnog_input", showWarnings = FALSE)

# -----------------------------
# 1. Read protein FASTA
# -----------------------------

faa_lines <- readLines("01_reference/complete.proteins.faa", warn = FALSE)

header_idx <- grep("^>", faa_lines)

if (length(header_idx) == 0) {
  stop("No FASTA headers found in complete.proteins.faa")
}

record_start <- header_idx
record_end <- c(header_idx[-1] - 1, length(faa_lines))

# -----------------------------
# 2. Parse FASTA records
# -----------------------------

record_list <- vector("list", length(record_start))

for (i in seq_along(record_start)) {
  
  header <- faa_lines[record_start[i]]
  
  if ((record_start[i] + 1) <= record_end[i]) {
    seq_lines <- faa_lines[(record_start[i] + 1):record_end[i]]
  } else {
    seq_lines <- character(0)
  }
  
  sequence <- paste(seq_lines, collapse = "")
  sequence <- gsub("\\s+", "", sequence)
  
  # First token after >
  first_token <- sub("^>", "", header)
  first_token <- sub(" .*", "", first_token)
  
  # Protein ID after final |
  # Example:
  # gnl|WGS:ZZZZ|egapxtmp_014493-P1 -> egapxtmp_014493-P1
  protein_id <- sub(".*\\|", "", first_token)
  
  # Gene ID without protein isoform suffix
  # Example:
  # egapxtmp_014493-P1 -> egapxtmp_014493
  gene_id <- sub("-P[0-9]+$", "", protein_id)
  
  record_list[[i]] <- data.frame(
    header = header,
    protein_id = protein_id,
    gene_id = gene_id,
    protein_length_aa = nchar(sequence),
    sequence = sequence,
    stringsAsFactors = FALSE
  )
}

records <- do.call(rbind, record_list)

# Remove empty or malformed gene IDs, just in case
records <- records[!is.na(records$gene_id) & records$gene_id != "", ]

# -----------------------------
# 3. Select longest isoform per gene
# -----------------------------

# Sort by gene ID, longest protein first, then protein_id for deterministic tie handling
records <- records[
  order(records$gene_id, -records$protein_length_aa, records$protein_id),
]

representative_records <- records[!duplicated(records$gene_id), ]

# Check genes with multiple isoforms
isoform_count_per_gene <- table(records$gene_id)
genes_with_multiple_isoforms <- sum(isoform_count_per_gene > 1)

# Check tied longest isoforms
max_len_by_gene <- ave(records$protein_length_aa, records$gene_id, FUN = max)
longest_candidates <- records[records$protein_length_aa == max_len_by_gene, ]
longest_candidate_count <- table(longest_candidates$gene_id)
genes_with_tied_longest_isoforms <- sum(longest_candidate_count > 1)

# -----------------------------
# 4. Safe FASTA wrapping function
# -----------------------------

wrap_sequence <- function(sequence, width = 60) {
  
  sequence_length <- nchar(sequence)
  
  if (is.na(sequence_length) || sequence_length == 0) {
    return(character(0))
  }
  
  starts <- seq(1, sequence_length, by = width)
  ends <- pmin(starts + width - 1, sequence_length)
  
  substring(sequence, starts, ends)
}

# -----------------------------
# 5. Write representative FASTA
# -----------------------------

output_faa <- "10_functional_annotation/eggnog_input/Qerythrinae_longest_protein_per_gene.faa"

fasta_out <- character()

for (i in seq_len(nrow(representative_records))) {
  
  gene_id <- representative_records$gene_id[i]
  protein_id <- representative_records$protein_id[i]
  sequence <- representative_records$sequence[i]
  protein_length <- representative_records$protein_length_aa[i]
  
  new_header <- paste0(
    ">",
    gene_id,
    " original_protein=",
    protein_id,
    " length_aa=",
    protein_length
  )
  
  seq_split <- wrap_sequence(sequence, width = 60)
  
  fasta_out <- c(fasta_out, new_header, seq_split)
}

writeLines(fasta_out, output_faa)

# -----------------------------
# 6. Save selection tables
# -----------------------------

# Save selected representative table without full sequences
representative_table <- representative_records[, c(
  "protein_id",
  "gene_id",
  "protein_length_aa",
  "header"
)]

write.table(
  representative_table,
  file = "10_functional_annotation/eggnog_input/longest_protein_per_gene_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Save all isoform lengths, useful for checking
all_isoform_table <- records[, c(
  "protein_id",
  "gene_id",
  "protein_length_aa",
  "header"
)]

write.table(
  all_isoform_table,
  file = "10_functional_annotation/eggnog_input/all_protein_isoform_lengths.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 7. Compare with protein-coding genes in GFF
# -----------------------------

gff_lines <- readLines("01_reference/complete.genomic.gff", warn = FALSE)

protein_coding_genes_in_gff <- sum(
  grepl("\tgene\t", gff_lines) &
    grepl("gene_biotype=protein_coding", gff_lines)
)

representative_fasta_header_count <- sum(grepl("^>", readLines(output_faa, warn = FALSE)))

# -----------------------------
# 8. Summary
# -----------------------------

representative_summary <- data.frame(
  item = c(
    "Total protein isoforms in complete.proteins.faa",
    "Unique genes parsed from protein FASTA",
    "Genes with more than one protein isoform",
    "Genes with tied longest isoforms",
    "Representative proteins written",
    "Representative FASTA header count",
    "Protein-coding genes in GFF",
    "Representative proteins minus protein-coding genes in GFF"
  ),
  count = c(
    nrow(records),
    length(unique(records$gene_id)),
    genes_with_multiple_isoforms,
    genes_with_tied_longest_isoforms,
    nrow(representative_records),
    representative_fasta_header_count,
    protein_coding_genes_in_gff,
    nrow(representative_records) - protein_coding_genes_in_gff
  )
)

write.table(
  representative_summary,
  file = "10_functional_annotation/eggnog_input/representative_protein_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 9. Print checks
# -----------------------------

representative_summary

head(representative_table, 10)

list.files("10_functional_annotation/eggnog_input")