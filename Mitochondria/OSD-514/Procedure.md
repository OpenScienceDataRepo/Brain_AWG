# Procedure Overview

To analyze NASA's OSD-514 Study, differential expression analysis was first performed using DESeq2 in R, yielding several differentially expressed genes. Subsequent GSEA revealed some enriched mitochondrial pathways, with only marginal evidence of mitochondrial involvement. It was concluded that while spaceflight induces detectable transcriptional changes, the dataset provides limited support for robust mitochondrial dysfunction at the pathway level.

# DESeq2 Analysis

```
# 0) Load libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(biomaRt)
})

# 1) Set working directory
setwd("~/Desktop/ADBR Mito/OSD-514/")  # change as needed

# 2) Load counts and metadata
counts <- read.csv("counts.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv", check.names = FALSE, header = TRUE)

# 3) Prepare metadata
rownames(metadata) <- trimws(metadata[,1])
metadata <- metadata[,-1]

# 4) Rename metadata columns
colnames(metadata)[colnames(metadata) == "Factor Value[Sex]"] <- "sex"
colnames(metadata)[colnames(metadata) == "Factor Value[Spaceflight]"] <- "condition"

# 5) Standardize counts column names
colnames(counts) <- toupper(trimws(colnames(counts)))
colnames(counts) <- sub("_CRRA.*", "", colnames(counts))
colnames(counts) <- gsub("_", " ", colnames(counts))

# 6) Standardize metadata rownames
rownames(metadata) <- toupper(gsub("_", " ", trimws(rownames(metadata))))

# 7) Map counts column IDs to metadata format
map_ids <- function(x){
  x <- gsub("^UG M", "SPACEFLIGHT MICROGRAVITY MALE ", x)
  x <- gsub("^UG F", "SPACEFLIGHT MICROGRAVITY FEMALE ", x)
  x <- gsub("^S 1G M", "SPACEFLIGHT 1G MALE ", x)
  x <- gsub("^S 1G F", "SPACEFLIGHT 1G FEMALE ", x)
  x <- gsub("^AGC M", "EARTH MALE ", x)
  x <- gsub("^AGC F", "EARTH FEMALE ", x)
  trimws(x)
}
colnames(counts) <- map_ids(colnames(counts))

# 8) Map metadata rownames IDs to the same format
map_metadata_ids <- function(x){
  x <- gsub(" M$", " MALE", x)
  x <- gsub(" F$", " FEMALE", x)
  trimws(x)
}
rownames(metadata) <- map_metadata_ids(rownames(metadata))

# 9) Check alignment
missing_in_metadata <- setdiff(colnames(counts), rownames(metadata))
missing_in_counts <- setdiff(rownames(metadata), colnames(counts))
if(length(missing_in_metadata) > 0 | length(missing_in_counts) > 0){
  stop("Mismatch between counts columns and metadata rownames. Inspect names!")
} else {
  message("Counts columns and metadata rownames match!")
}

# 10) Reorder counts columns to match metadata
counts <- counts[, rownames(metadata)]

# 11) Convert counts to integers
counts <- round(counts)
counts <- apply(counts, 2, as.integer)
counts <- as.data.frame(counts)
rownames(counts) <- rownames(read.csv("counts.csv", row.names = 1, check.names = FALSE))

# 12) Optional: annotate genes with biomaRt (Drosophila)
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(counts),
  mart = ensembl
)
counts$gene_symbol <- genes$external_gene_name[match(rownames(counts), genes$ensembl_gene_id)]
counts$gene_symbol[is.na(counts$gene_symbol)] <- rownames(counts)[is.na(counts$gene_symbol)]
rownames(counts) <- make.unique(counts$gene_symbol)
counts$gene_symbol <- NULL

# 13) Collapse conditions into main groups
metadata$condition_group <- rownames(metadata) # use rownames as source
metadata$condition_group <- gsub("^SPACEFLIGHT MICROGRAVITY.*", "SPACEFLIGHT_MICROGRAVITY", metadata$condition_group)
metadata$condition_group <- gsub("^SPACEFLIGHT 1G.*", "SPACEFLIGHT_1G", metadata$condition_group)
metadata$condition_group <- gsub("^EARTH.*", "EARTH", metadata$condition_group)

# 14) Make factor with EARTH as reference
metadata$condition_group <- factor(metadata$condition_group,
                                   levels = c("EARTH", "SPACEFLIGHT_1G", "SPACEFLIGHT_MICROGRAVITY"))

# 15) Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata[, c("sex", "condition_group")],
  design = ~ sex + condition_group
)

# 16) Prefilter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# 17) Run DESeq2
dds <- DESeq(dds)

# 18) Check available results names
resultsNames(dds)
# Should show:
# "Intercept" "sex_MALE_vs_FEMALE" "condition_group_SPACEFLIGHT_1G_vs_EARTH" "condition_group_SPACEFLIGHT_MICROGRAVITY_vs_EARTH"

# 19) Extract results using 'coef'
# Microgravity vs Earth
res_MG_vs_E <- lfcShrink(dds,
                         coef="condition_group_SPACEFLIGHT_MICROGRAVITY_vs_EARTH",
                         type="apeglm")
write.csv(as.data.frame(res_MG_vs_E), "DESeq2_results_Microgravity_vs_Earth.csv")

# 1G vs Earth
res_1G_vs_E <- lfcShrink(dds,
                         coef="condition_group_SPACEFLIGHT_1G_vs_EARTH",
                         type="apeglm")
write.csv(as.data.frame(res_1G_vs_E), "DESeq2_results_1G_vs_Earth.csv")

```
# Volcano Plot Generation
```

# Function to create volcano plot
make_volcano <- function(res, title, lfc_cutoff = 1, p_cutoff = 0.05) {
  # Drop rows with NA values
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

  # Add significance category
  res$significance <- "Not significant"
  res$significance[res$padj < p_cutoff & res$log2FoldChange > lfc_cutoff] <- "Upregulated"
  res$significance[res$padj < p_cutoff & res$log2FoldChange < -lfc_cutoff] <- "Downregulated"

  # Compute -log10(padj)
  res$negLogP <- -log10(res$padj)

  # Volcano plot (return object)
  p <- ggplot(res, aes(x = log2FoldChange, y = negLogP, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Upregulated" = "red",
                                  "Downregulated" = "blue",
                                  "Not significant" = "grey")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
    labs(title = title,
         x = "log2(Fold Change)",
         y = "-log10(Adjusted p-value)") +
    theme_minimal()

  return(p)
}

# --- Generate and save plots ---
p1 <- make_volcano(res_MG_vs_E, "Microgravity vs Earth")
ggsave("Volcano_Microgravity_vs_Earth.png", p1, width = 7, height = 6, dpi = 300)

p2 <- make_volcano(res_1G_vs_E, "1G vs Earth")
ggsave("Volcano_1G_vs_Earth.png", p2, width = 7, height = 6, dpi = 300)

```
# Pre-GSEA File Creation & Preparation
A ranked gene list for GSEA was created using the below code
```
# 0) Set working directory
setwd("~/Desktop/ADBR Mito/OSD-514/RESULTS")  # change as needed

# 1) Load required libraries
library(dplyr)

# 2) Function to generate .rnk file
generate_rnk <- function(file_in, file_out) {
  # Read CSV with rownames
  res <- read.csv(file_in, header = TRUE, check.names = FALSE, row.names = 1)

  # Add rownames as a proper column
  res$gene <- rownames(res)

  # Filter out rows with NA in log2FoldChange or gene
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$gene) & res$gene != "", ]

  # Use log2FoldChange for ranking
  ranked_genes <- res$log2FoldChange
  names(ranked_genes) <- res$gene

  # Order decreasing
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  # Save as .rnk for GSEA
  write.table(
    data.frame(Gene = names(ranked_genes), Score = ranked_genes),
    file_out,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  message(paste("Saved .rnk file:", file_out))
}

# 3) Apply function to comparisons

# Microgravity vs Earth
generate_rnk(
  "DESeq2_results_Microgravity_vs_Earth.csv",
  "DEG_ranked_list_Microgravity_vs_Earth.rnk"
)

# 1G in centrifuge vs Earth
generate_rnk(
  "DESeq2_results_1G_vs_Earth.csv",
  "DEG_ranked_list_1G_vs_Earth.rnk"
)
```
I also used Jack Cheng's gmt file (from GeneMatrix on GitHub) for GSEA

Initial GSEA yielded an error, likely a formatting error detected by the below code:
```
# 0 Load libraries
library(dplyr)
library(stringr)


# 1 Set working directory
setwd("~/Desktop/ADBR Mito/OSD-514/RESULTS")

# 2 Load your ranked gene list (from DESeq2 or similar)
# Must have at least two columns: "Gene" and "Stat" (log2FC, t-statistic, etc.)
ranked_file <- "DEG_ranked_list.rnk" # replace with your file

# Use base R to read it (no readr needed)
ranked_df <- read.delim(ranked_file, header = TRUE, stringsAsFactors = FALSE)

# Example columns: Gene, log2FoldChange

# 3 Load GMT file
gmt_file <- "Drosophila_KEGG.gmt" # replace with actual GMT filename

# Read GMT: each line is one gene set
gmt_lines <- readLines(gmt_file)
gmt_list <- lapply(gmt_lines, function(line) {
  parts <- str_split(line, "\t")[[1]]
  # First element = GeneSet name, 2nd = description, rest = genes
  list(name = parts[1], desc = parts[2], genes = parts[-c(1,2)])
})

# 4 Check gene ID overlap
overlap_results <- data.frame(
  GeneSet = sapply(gmt_list, function(x) x$name),
  SetSize = sapply(gmt_list, function(x) length(x$genes)),
  Overlap = sapply(gmt_list, function(x) sum(x$genes %in% ranked_df$Gene))
)

# 5 Report summary
cat("Total gene sets:", nrow(overlap_results), "\n")
cat("Gene sets with ≥1 overlapping gene:", sum(overlap_results$Overlap >= 1), "\n")
cat("Gene sets passing typical min size (>=15) threshold:",
    sum(overlap_results$Overlap >= 15), "\n")

# Optional: view top overlaps
top_overlaps <- overlap_results %>% arrange(desc(Overlap)) %>% head(20)
print(top_overlaps)

# 6 Identify problematic gene sets
problematic <- overlap_results[overlap_results$Overlap == 0, ]
if(nrow(problematic) > 0){
  cat("Warning: The following gene sets have ZERO overlap with your ranked list:\n")
  print(problematic$GeneSet)
} else {
  cat("All gene sets have at least one overlapping gene.\n")
}
```
The below code was used to correct and replace the rnk file
```
# 0 Load libraries
library(dplyr)
library(stringr)
library(biomaRt)

# 1 Set working directory
setwd("~/Desktop/ADBR Mito/OSD-514/RESULTS")  # adjust as needed

# 2 Load your ranked gene list
# Replace with your file name
ranked_file <- "DEG_ranked_list.rnk"

# Load as plain table; set column names
ranked_df <- read.delim(ranked_file, header = FALSE, stringsAsFactors = FALSE)
colnames(ranked_df) <- c("Gene", "Stat")

# Quick check
head(ranked_df)

# 3 Map FlyBase IDs to gene symbols (for KEGG GMT)
fly <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

mapping <- getBM(
  attributes = c("flybase_gene_id", "external_gene_name"),
  filters = "flybase_gene_id",
  values = ranked_df$Gene,
  mart = fly
)

# Merge mapped symbols back into ranked list
ranked_df <- merge(ranked_df, mapping, by.x = "Gene", by.y = "flybase_gene_id")
ranked_df$Gene <- ranked_df$external_gene_name

# Remove any rows without mapped symbols
ranked_df <- ranked_df[!is.na(ranked_df$Gene), ]

head(ranked_df)

# 4 Load GMT file
gmt_file <- "Drosophila_KEGG.gmt"  # replace with actual GMT filename
gmt_lines <- readLines(gmt_file)

gmt_list <- lapply(gmt_lines, function(line) {
  parts <- str_split(line, "\t")[[1]]
  list(
    name = parts[1],
    desc = parts[2],
    genes = parts[-c(1,2)]
  )
})

# 5 Check gene ID overlap
overlap_results <- data.frame(
  GeneSet = sapply(gmt_list, function(x) x$name),
  SetSize = sapply(gmt_list, function(x) length(x$genes)),
  Overlap = sapply(gmt_list, function(x) sum(x$genes %in% ranked_df$Gene))
)

# 6 Report summary
cat("Total gene sets:", nrow(overlap_results), "\n")
cat("Gene sets with ≥1 overlapping gene:", sum(overlap_results$Overlap >= 1), "\n")
cat("Gene sets passing typical min size (>=15) threshold:", sum(overlap_results$Overlap >= 15), "\n")

# Optional: view top overlaps
top_overlaps <- overlap_results %>% arrange(desc(Overlap)) %>% head(20)
print(top_overlaps)

# Identify problematic gene sets
problematic <- overlap_results[overlap_results$Overlap == 0, ]
if(nrow(problematic) > 0){
  cat("Warning: The following gene sets have ZERO overlap with your ranked list:\n")
  print(problematic$GeneSet)
} else {
  cat("All gene sets have at least one overlapping gene.\n")
}

write.table(ranked_df[, c("Gene", "Stat")],
            file = "DEG_ranked_symbols.rnk",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# Save corrected rnk file
write.table(ranked_df[, c("Gene", "Stat")],
            file = "DEG_ranked_symbols.rnk",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

```
# GSEA Analysis and Results
GSEA was run on both conditions (1G vs Earth, Microgravity vs Earth) in order to investigate more subtle changes and overall trends in gene expression.

The Enrichment results were saved and analyzed (available in the RESULTS folder)

# Heatmap Generation
Using OSD-514's normalized counts list, heatmaps were generated for both comparisons using the below code.

NOTE: a validated id list from flybase was created to harmonize the databases

```
# 1) Set working directory
setwd("~/Desktop/ADBR Mito/JPT TEST/DESeq2")  # adjust as needed

# 2) Load required packages
packages <- c("pheatmap", "RColorBrewer", "dplyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# 3) Load normalized counts and DESeq2 results
norm_counts <- read.csv("norm_counts.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
res_MG_vs_E <- read.csv("DESeq2_results_Microgravity_vs_Earth.csv", row.names = 1, stringsAsFactors = FALSE)
res_1G_vs_E <- read.csv("DESeq2_results_1G_vs_Earth.csv", row.names = 1, stringsAsFactors = FALSE)

  # Ensure counts rownames are character
rownames(norm_counts) <- trimws(as.character(rownames(norm_counts)))

# 4) Load FlyBase mapping table (FBgn → Symbol)
gene_map <- read.delim("flybase_mapping.csv", header = FALSE, stringsAsFactors = FALSE)

  # Fix mapping table (assume FBgn is in col2, Symbol in col3)
colnames(gene_map) <- c("col1", "FBgn", "Symbol")
gene_map <- gene_map[, c("FBgn", "Symbol")]
gene_map$FBgn <- trimws(as.character(gene_map$FBgn))
gene_map$Symbol <- trimws(as.character(gene_map$Symbol))

  # Remove first row if it contains headers
if(gene_map$FBgn[1] == "validated_id") {
  gene_map <- gene_map[-1, ]
}

# 5) Map counts rownames to symbols
mapped_symbols <- gene_map$Symbol[match(rownames(norm_counts), gene_map$FBgn)]
valid_idx <- !is.na(mapped_symbols)
norm_counts <- norm_counts[valid_idx, , drop = FALSE]
rownames(norm_counts) <- mapped_symbols[valid_idx]

cat("Number of genes successfully mapped:", nrow(norm_counts), "\n")
cat("Mapped counts rownames (first 10):\n")
print(head(rownames(norm_counts), 10))

# 6) Heatmap function with max_display & readable labels
generate_heatmap <- function(norm_counts, res_obj, top_n = 30, max_display = 50,
                             main_title = "Heatmap", file_name = "heatmap.jpg") {

  # Filter significant genes and sort by padj
  top_genes <- res_obj %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj)

  # Limit to top_n
  top_genes <- head(top_genes, top_n)
  top_gene_names <- rownames(top_genes)
  cat("Number of genes considered before filtering:", length(top_gene_names), "\n")

  # Subset counts
  genes_exist <- intersect(top_gene_names, rownames(norm_counts))
  if(length(genes_exist) == 0) stop("No top genes exist in normalized counts. Check your mapping!")

  max_genes <- min(length(genes_exist), max_display)
  heatmap_data <- norm_counts[genes_exist[1:max_genes], , drop = FALSE]

  # Remove rows with NA/NaN/Inf or zero variance
  heatmap_data <- heatmap_data[apply(heatmap_data, 1, function(x) all(is.finite(x)) && sd(x) != 0), , drop = FALSE]
  if(nrow(heatmap_data) < 2) stop("Not enough valid genes to generate heatmap.")

  # Scale rows
  heatmap_data <- t(scale(t(heatmap_data)))

  # Shorten sample names for readability
  colnames(heatmap_data) <- sub("_CRRA.*", "", colnames(heatmap_data))

  # Colors
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)

  # Generate and save heatmap
  jpeg(file_name, width = 2000, height = 1400, res = 200)
  pheatmap(heatmap_data,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           color = colors,
           main = main_title)
  dev.off()

  cat("Heatmap saved:", file_name, "\n")
  cat("Number of genes plotted:", nrow(heatmap_data), "\n")
}

# 7) Function to select samples by prefixes
subset_samples_by_prefix <- function(norm_counts, prefixes) {
  cols <- colnames(norm_counts)
  matched <- unlist(lapply(prefixes, function(p) grep(p, cols, value = TRUE)))
  if(length(matched) < 2) stop("Less than 2 samples matched the patterns.")
  return(matched)
}

# 8) Subset samples for comparisons
samples_MG_E <- subset_samples_by_prefix(norm_counts, c("ug_", "AGC_"))
norm_counts_MG_E <- norm_counts[, samples_MG_E, drop = FALSE]

samples_1G_E <- subset_samples_by_prefix(norm_counts, c("s_1g_", "AGC_"))
norm_counts_1G_E <- norm_counts[, samples_1G_E, drop = FALSE]

# 9) Optional: check number of significant genes
cat("Significant genes (padj < 0.05) Microgravity vs Earth:",
    sum(!is.na(res_MG_vs_E$padj) & res_MG_vs_E$padj < 0.05), "\n")
cat("Significant genes (padj < 0.05) 1G vs Earth:",
    sum(!is.na(res_1G_vs_E$padj) & res_1G_vs_E$padj < 0.05), "\n")

# 10) Generate heatmaps
generate_heatmap(norm_counts_MG_E, res_MG_vs_E, top_n = 30, max_display = 30,
                 main_title = "Top DE Genes: Microgravity vs Earth",
                 file_name = "Heatmap_Microgravity_vs_Earth.jpg")

generate_heatmap(norm_counts_1G_E, res_1G_vs_E, top_n = 30, max_display = 30,
                 main_title = "Top DE Genes: 1G vs Earth",
                 file_name = "Heatmap_1G_vs_Earth.jpg")
```
