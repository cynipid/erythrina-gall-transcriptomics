# 01_setup_check.R
# Purpose:
# Check whether R can correctly read the gene count matrix and sample metadata.

# Set project folder
setwd("//wsl.localhost/Ubuntu/home/dodi/projects/qe_genome_guided")

# Read gene count matrix
counts <- read.delim("06_counts/gene_counts_matrix.tsv",
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

# Read sample metadata
metadata <- read.delim("06_counts/sample_metadata.tsv",
                       header = TRUE,
                       row.names = 1,
                       check.names = FALSE)

# Check dimensions
dim(counts)
dim(metadata)

# Check sample names
colnames(counts)
rownames(metadata)

# Check whether count matrix columns match metadata rows
all(colnames(counts) == rownames(metadata))

# Look at first few rows
head(counts)

# Show metadata
metadata