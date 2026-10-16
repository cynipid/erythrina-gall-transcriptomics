# 10_KEGG_pathway_enrichment_primary.R
# KEGG pathway over-representation analysis
# for official insect DEG sets
#
# DEG definition:
# padj < 0.01 and |log2FC| >= 1
#
# KEGG enrichment significance:
# BH-adjusted P < 0.05
#
# Only canonical mapXXXXX pathways are used.

setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")


# ------------------------------------------------------------
# 1. Create output folder
# ------------------------------------------------------------

dir.create(
  "10_functional_annotation/enrichment_results_primary",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Read complete DESeq2 results
# ------------------------------------------------------------

all_res <- read.delim(
  "08_results_primary/deseq2_pupa_vs_larva_all_genes.tsv",
  header = TRUE,
  check.names = FALSE
)

# Background = DESeq2-tested genes with valid adjusted P-value

all_res_valid <- all_res[
  !is.na(all_res$padj),
]

background_genes <- unique(
  all_res_valid$Geneid
)


# ------------------------------------------------------------
# 3. Read official DEG lists
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
# 4. Read KEGG pathway mapping
# ------------------------------------------------------------

gene_to_KEGG <- read.delim(
  "10_functional_annotation/parsed_mappings/gene_to_KEGG_pathway.tsv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Keep canonical map pathways only
# ------------------------------------------------------------

gene_to_KEGG_map <- gene_to_KEGG[
  grepl("^map", gene_to_KEGG$KEGG_Pathway),
]

gene_to_KEGG_map <- unique(
  gene_to_KEGG_map[, c("Geneid", "KEGG_Pathway")]
)
# ------------------------------------------------------------
# Keep only pathways present in the current KEGG database
# ------------------------------------------------------------

library(KEGGREST)

current_KEGG <- KEGGREST::keggList("pathway")

current_KEGG_ids <- sub(
  "^path:",
  "",
  names(current_KEGG)
)

gene_to_KEGG_map <- gene_to_KEGG_map[
  gene_to_KEGG_map$KEGG_Pathway %in% current_KEGG_ids,
]

# ------------------------------------------------------------
# 6. KEGG ORA function
# ------------------------------------------------------------

run_KEGG_ORA <- function(
    gene_list,
    background,
    term_mapping,
    min_term_size = 3,
    max_term_size = 500
) {
  
  # Restrict mapping to DESeq2 background
  
  mapping_bg <- term_mapping[
    term_mapping$Geneid %in% background,
  ]
  
  mapping_bg <- mapping_bg[
    , c("Geneid", "KEGG_Pathway")
  ]
  
  colnames(mapping_bg) <- c(
    "Geneid",
    "term"
  )
  
  mapping_bg <- unique(mapping_bg)
  
  
  # Restrict DEG list to DESeq2 background
  
  gene_list <- intersect(
    unique(gene_list),
    background
  )
  
  
  # KEGG-annotated background genes
  
  annotated_background <- unique(
    mapping_bg$Geneid
  )
  
  
  # DEGs with KEGG pathway annotation
  
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
  
  
  # Split genes by pathway
  
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
      
      
      # Keep pathways represented by 3-500 genes
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
      # Zero-overlap pathways are retained (p = 1)
      
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
  
  
  # Benjamini-Hochberg FDR correction
  # across all eligible KEGG pathways
  
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
# 7. Run KEGG pathway enrichment
# ------------------------------------------------------------

KEGG_pupa <- run_KEGG_ORA(
  gene_list = pupa_up_genes,
  background = background_genes,
  term_mapping = gene_to_KEGG_map,
  min_term_size = 3,
  max_term_size = 500
)

KEGG_larva <- run_KEGG_ORA(
  gene_list = larva_up_genes,
  background = background_genes,
  term_mapping = gene_to_KEGG_map,
  min_term_size = 3,
  max_term_size = 500
)


# ------------------------------------------------------------
# 8. Significant pathways
# ------------------------------------------------------------

KEGG_pupa_sig <- KEGG_pupa[
  KEGG_pupa$padj < 0.05,
]

KEGG_larva_sig <- KEGG_larva[
  KEGG_larva$padj < 0.05,
]


# ------------------------------------------------------------
# 9. Save complete pathway tables
# ------------------------------------------------------------

write.table(
  KEGG_pupa,
  "10_functional_annotation/enrichment_results_primary/KEGG_pathway_enrichment_pupa_up.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  KEGG_larva,
  "10_functional_annotation/enrichment_results_primary/KEGG_pathway_enrichment_larva_up.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Save significant pathway tables
# ------------------------------------------------------------

write.table(
  KEGG_pupa_sig,
  "10_functional_annotation/enrichment_results_primary/KEGG_pathway_enrichment_pupa_up_FDR005.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  KEGG_larva_sig,
  "10_functional_annotation/enrichment_results_primary/KEGG_pathway_enrichment_larva_up_FDR005.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Short summary
# ------------------------------------------------------------

cat("\nKEGG pathway enrichment\n")

cat("\nPupa-up\n")
cat("DEGs:", length(pupa_up_genes), "\n")
cat("DEGs with KEGG pathway:",
    KEGG_pupa$gene_list_annotated[1], "\n")
cat("Pathways tested:",
    nrow(KEGG_pupa), "\n")
cat("Significant pathways:",
    nrow(KEGG_pupa_sig), "\n")

cat("\nLarva-up\n")
cat("DEGs:", length(larva_up_genes), "\n")
cat("DEGs with KEGG pathway:",
    KEGG_larva$gene_list_annotated[1], "\n")
cat("Pathways tested:",
    nrow(KEGG_larva), "\n")
cat("Significant pathways:",
    nrow(KEGG_larva_sig), "\n")