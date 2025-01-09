#!/usr/bin/env Rscript

# Load libraries
library(DESeq2)

# Input paths
count_data_file <- data/output/pair_end/test/counts/counts_summary.txt
metadata_file <- "data/metadata.tsv"

# Output paths
deseq2_results <- "data/output/pair_end/test/deseq2_results.csv"

# Load data
count_data <- read.table(count_data_file,header = FALSE, stringsAsFactors = FALSE)
metadata <- read.table(metadata_file, header=TRUE, sep="\t" )

# Clean data
cleaned_counts <- raw_counts[, -c(2:6)]
cleaned_counts <- cleaned_counts[-1, ]

# Assign new columns (samples)
colnames(cleaned_counts) <- c("Geneid", "ERR031027", "ERR031028", "ERR031029", "ERR031030", "ERR031031", "ERR031032", "ERR031033", "ERR299297")

# Assigning first column of cleaned_counts(manipulation to match tables)
rownames(cleaned_counts) <- cleaned_counts[, 1]
cleaned_counts <- cleaned_counts[, -1]

# Check if the column names in counts data match the 'sample' column in metadata
all(colnames(cleaned_counts) %in% metadata$sample)

# Rearrange the metadata rows to match the order of columns in the cleaned_counts table
metadata <- metadata[match(colnames(cleaned_counts), metadata$sample), ]

#Convert the 'group' column in metadata to a factor
# R will treat these as categorical data, rather than just character strings
# Tumor and non Tumor in this case
metadata$group <- factor(metadata$group)

# Convert data in the table into numeric
cleaned_counts[] <- lapply(cleaned_counts, as.numeric)

# Perform DESeq2
# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = cleaned_counts,
                              colData = metadata,
                              design = ~ group)

# Run the DESeq2 analysis
dds <- DESeq(dds)

# Get the results for the differential expression (e.g., for the 'group' variable)
res <- results(dds)

# Save results
write.csv(as.data.frame(res), file = deseq2_results)


