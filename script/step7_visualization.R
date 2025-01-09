#!/usr/bin/env Rscript

# Load libraries
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(biomaRt)
library(org.Hs.eg.db)

# Input path
deseq2_results <- "data/output/pair_end/test/deseq2_results.csv

# Output path
volcano_plot <- "data/output/pair_end/test/plots/volcano_plot.pdf"
heatmap_plot <- "data/output/pair_end/test/plots/heatmap.pdf"

# Load data
res <- read.csv(deseq2_results, row.names = 1)

# Remove rows with missing values in the relevant columns (`log2FoldChange` and `pvalue`)
res_clean <- res[complete.cases(res$log2FoldChange, res$pvalue), ]

# Filter the data to only include significantly different genes
filtered_res <- res_clean[abs(res_clean$log2FoldChange) > 1 & res_clean$pvalue < 0.05, ]

# Select top 40 genes based on the ranking
top_genes <- filtered_res[filtered_res$rank <= 40, ]

# Volcano Plot

# Define thresholds
x_threshold <- 1  # Example threshold for Log2 Fold Change (positive threshold)
x_negative_threshold <- -1  # Example threshold for Log2 Fold Change (negative threshold)
y_threshold <- -log10(0.05)  # Example threshold for p-value (adjusted)
y_lower_threshold <- 0  # Lower y-threshold for -log10(p-value)

# Assuming 'res_clean' is your data frame and it contains Ensembl gene IDs or other identifiers
# Ensure gene names are included. If not, retrieve them using AnnotationDbi
# For example, using org.Hs.eg.db to get external gene names from Ensembl IDs
res_clean$GeneName <- mapIds(org.Hs.eg.db, keys = rownames(res_clean), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Create a new column for colors based on your conditions
res_clean$color <- "grey"  # Default to grey
res_clean$color[res_clean$log2FoldChange > x_threshold & -log10(res_clean$pvalue) > y_threshold] <- "red"  # red for significant up-regulated
res_clean$color[res_clean$log2FoldChange <= x_negative_threshold & -log10(res_clean$pvalue) > y_threshold] <- "blue"  # blue for significant down-regulated
res_clean$color[res_clean$log2FoldChange > x_threshold & -log10(res_clean$pvalue) <= y_lower_threshold] <- "yellow"  # yellow for low p-value but no significant fold change
res_clean$color[res_clean$log2FoldChange <= x_negative_threshold & -log10(res_clean$pvalue) <= y_lower_threshold] <- "grey"  # grey for non-significant

# Sort by log2FoldChange and select top 20 genes
top_genes_upregulated <- res_clean[order(res_clean$log2FoldChange, decreasing = TRUE), ][1:40, ]
top_genes_downregulated <- res_clean[order(res_clean$log2FoldChange, decreasing = FALSE), ][1:40, ]

# Combine the top 20 up-regulated and top 20 down-regulated genes
top_genes <- rbind(top_genes_upregulated, top_genes_downregulated)

# Create the volcano plot
volcano_plot <- ggplot(res_clean, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = color), alpha = 0.5) +  # 50% transparency for dots
  scale_color_manual(
    values = c("blue" = "blue", "grey" = "grey", "red" = "red", "yellow" = "yellow"),
    labels = c("Down-regulated", "Non-significant", "Up-regulated", "Low p-value, no significant fold change")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank(),  # Remove legend title
        legend.position = "right") +  # Position the legend to the right
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 p-value") +
  # Add dotted threshold lines for y-axis (upper and lower) and x-axis (positive and negative)
  geom_hline(yintercept = y_threshold, linetype = "dotted", color = "black", size = 0.5) +  # Upper y-threshold
  geom_hline(yintercept = y_lower_threshold, linetype = "dotted", color = "black", size = 0.5) +  # Lower y-threshold
  geom_vline(xintercept = x_threshold, linetype = "dotted", color = "black", size = 0.5) +  # Positive x-threshold
  geom_vline(xintercept = x_negative_threshold, linetype = "dotted", color = "black", size = 0.5) +  # Negative x-threshold
  # Set custom axis limits for x and y axes
  scale_x_continuous(limits = c(min(res_clean$log2FoldChange, na.rm = TRUE), max(res_clean$log2FoldChange, na.rm = TRUE))) + 
  scale_y_continuous(limits = c(0, max(-log10(res_clean$pvalue), na.rm = TRUE) + 1)) + # Automatically adjust y-axis limit based on data
  # Label top 20 genes (both up-regulated and down-regulated) with small font size to avoid overlap
  geom_text_repel(aes(label = ifelse(GeneName %in% top_genes$GeneName, GeneName, "")), 
                  box.padding = 0.35, point.padding = 0.5,
                  segment.color = "grey50", max.overlaps = 500, 
                  size = 2) +  # Reduce font size to avoid overlap
  # Add annotation about the annotation database in the top-right corner
  annotate("text", x = Inf, y = Inf, label = "Annotation: org.Hs.eg.db", hjust = 1, vjust = 1, size = 3, angle = 0, color = "black")

# Save the volcano plot as PDF
ggsave(volcano_plot)


expression_data <- counts(dds, normalized = TRUE)

top_40_genes <- res_clean[order(abs(res_clean$log2FoldChange), decreasing = TRUE), ][1:40, ]

top_40_genes$GeneName <- mapIds(org.Hs.eg.db, 
                                keys = rownames(top_40_genes), 
                                column = "SYMBOL", 
                                keytype = "ENSEMBL", 
                                multiVals = "first")

colnames(expression_data) <- metadata$sample

annotation <- data.frame(Group = metadata$group)

rownames(annotation) <- metadata$sample

top_genes <- rownames(top_40_genes)

expression_data <- expression_data[top_genes, , drop = FALSE]

# Replace row names in expression_data with gene symbols from top_40_genes
rownames(expression_data) <- top_40_genes$GeneName

# Plot the heatmap
heatmap_plot <- pheatmap(
  expression_data,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation,
  show_rownames = TRUE,  # Display row names (gene names)
  show_colnames = TRUE,  # Display column names (sample names)
  main = "Heatmap of Top 40 Genes",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)

ggsave(heatmap_plot)



