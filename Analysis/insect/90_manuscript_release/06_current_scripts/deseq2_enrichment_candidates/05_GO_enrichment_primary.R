# 05_GO_enrichment_primary.R
# GO over-representation analysis for official insect DEGs
#
# DEG definition:
# padj < 0.01 and |log2FC| >= 1
#
# GO enrichment significance:
# BH-adjusted P < 0.05
#
# Only currently valid GO terms in GO.db are used.

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 1. Load GO database
# ------------------------------------------------------------

library(GO.db)
library(AnnotationDbi)


# ------------------------------------------------------------
# 2. Create output folder
# ------------------------------------------------------------

dir.create(
  "10_functional_annotation/enrichment_results_primary",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Read complete DESeq2 results
# ------------------------------------------------------------

all_res <- read.delim(
  "08_results_primary/deseq2_pupa_vs_larva_all_genes.tsv",
  header = TRUE,
  check.names = FALSE
)

# Background = genes that received a valid DESeq2 adjusted P-value

all_res_valid <- all_res[
  !is.na(all_res$padj),
]

background_genes <- unique(
  all_res_valid$Geneid
)


# ------------------------------------------------------------
# 4. Read official DEG lists
# ------------------------------------------------------------

pupa_up <- read.delim(
  "08_results_primary/DEGs_pupa_up.tsv",
  header = TRUE,
  check.names = FALSE
)

larva_up <- read.delim(
  "08_results_primary/DEGs_larva_up.tsv",
  header = TRUE,
  check.names = FALSE
)

pupa_up_genes <- unique(
  pupa_up$Geneid
)

larva_up_genes <- unique(
  larva_up$Geneid
)


# ------------------------------------------------------------
# 5. Read GO mapping
# ------------------------------------------------------------

gene_to_GO <- read.delim(
  "10_functional_annotation/parsed_mappings/gene_to_GO.tsv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 6. Keep only current GO terms
# ------------------------------------------------------------

valid_GO_ids <- AnnotationDbi::keys(
  GO.db,
  keytype = "GOID"
)

gene_to_GO_current <- gene_to_GO[
  gene_to_GO$GO %in% valid_GO_ids,
]

gene_to_GO_current <- unique(
  gene_to_GO_current[, c("Geneid", "GO")]
)


# ------------------------------------------------------------
# 7. GO over-representation function
# ------------------------------------------------------------

run_GO_ORA <- function(
    gene_list,
    background,
    term_mapping,
    min_term_size = 5,
    max_term_size = 500
) {
  
  # Mapping restricted to tested DESeq2 background
  mapping_bg <- term_mapping[
    term_mapping$Geneid %in% background,
  ]
  
  mapping_bg <- mapping_bg[
    , c("Geneid", "GO")
  ]
  
  colnames(mapping_bg) <- c(
    "Geneid",
    "term"
  )
  
  mapping_bg <- unique(mapping_bg)
  
  # DEG list restricted to tested background
  gene_list <- intersect(
    unique(gene_list),
    background
  )
  
  # Genes with valid GO annotation
  annotated_background <- unique(
    mapping_bg$Geneid
  )
  
  gene_list_annotated <- intersect(
    gene_list,
    annotated_background
  )
  
  N <- length(
    annotated_background
  )
  
  n <- length(
    gene_list_annotated
  )
  
  # Split genes by GO term
  term_list <- split(
    mapping_bg$Geneid,
    mapping_bg$term
  )
  
  results <- lapply(
    names(term_list),
    function(term) {
      
      term_genes <- unique(
        term_list[[term]]
      )
      
      K <- length(
        term_genes
      )
      
      # Keep GO terms represented by 5-500 genes
      # in the annotated background
      
      if (
        K < min_term_size ||
        K > max_term_size
      ) {
        return(NULL)
      }
      
      overlapping_genes <- intersect(
        gene_list_annotated,
        term_genes
      )
      
      x <- length(
        overlapping_genes
      )
      
      # Hypergeometric over-representation test
      # Terms with zero overlap are retained and receive p = 1
      
      pvalue <- phyper(
        q = x - 1,
        m = K,
        n = N - K,
        k = n,
        lower.tail = FALSE
      )
      
      data.frame(
        term = term,
        overlap = x,
        gene_list_annotated = n,
        term_size_in_background = K,
        annotated_background_size = N,
        pvalue = pvalue,
        genes = paste(
          overlapping_genes,
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- do.call(
    rbind,
    results
  )
  
  if (
    is.null(results) ||
    nrow(results) == 0
  ) {
    return(data.frame())
  }
  
  # Benjamini-Hochberg correction across
  # all eligible GO terms
  
  results$padj <- p.adjust(
    results$pvalue,
    method = "BH"
  )
  
  results$GeneRatio <- paste0(
    results$overlap,
    "/",
    results$gene_list_annotated
  )
  
  results$BgRatio <- paste0(
    results$term_size_in_background,
    "/",
    results$annotated_background_size
  )
  
  results <- results[
    order(
      results$padj,
      results$pvalue
    ),
  ]
  
  return(results)
}


# ------------------------------------------------------------
# 8. Run GO enrichment
# ------------------------------------------------------------

GO_pupa <- run_GO_ORA(
  gene_list = pupa_up_genes,
  background = background_genes,
  term_mapping = gene_to_GO_current,
  min_term_size = 5,
  max_term_size = 500
)

GO_larva <- run_GO_ORA(
  gene_list = larva_up_genes,
  background = background_genes,
  term_mapping = gene_to_GO_current,
  min_term_size = 5,
  max_term_size = 500
)


# ------------------------------------------------------------
# 9. Significant GO terms
# ------------------------------------------------------------

GO_pupa_sig <- GO_pupa[
  GO_pupa$padj < 0.05,
]

GO_larva_sig <- GO_larva[
  GO_larva$padj < 0.05,
]


# ------------------------------------------------------------
# 10. Save complete GO tables
# ------------------------------------------------------------

write.table(
  GO_pupa,
  "10_functional_annotation/enrichment_results_primary/GO_enrichment_pupa_up.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  GO_larva,
  "10_functional_annotation/enrichment_results_primary/GO_enrichment_larva_up.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Save significant GO tables
# ------------------------------------------------------------

write.table(
  GO_pupa_sig,
  "10_functional_annotation/enrichment_results_primary/GO_enrichment_pupa_up_FDR005.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  GO_larva_sig,
  "10_functional_annotation/enrichment_results_primary/GO_enrichment_larva_up_FDR005.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Short summary
# ------------------------------------------------------------

cat("\nCurrent GO analysis\n")

cat("Original GO mappings:",
    nrow(gene_to_GO), "\n")

cat("Current GO mappings:",
    nrow(gene_to_GO_current), "\n")

cat("\nPupa-up\n")
cat("DEGs:", length(pupa_up_genes), "\n")
cat("DEGs with current GO:",
    GO_pupa$gene_list_annotated[1], "\n")
cat("GO terms tested:",
    nrow(GO_pupa), "\n")
cat("Significant GO terms:",
    nrow(GO_pupa_sig), "\n")

cat("\nLarva-up\n")
cat("DEGs:", length(larva_up_genes), "\n")
cat("DEGs with current GO:",
    GO_larva$gene_list_annotated[1], "\n")
cat("GO terms tested:",
    nrow(GO_larva), "\n")
cat("Significant GO terms:",
    nrow(GO_larva_sig), "\n")