# ============================================================
# Emory Drosophila Tau Model - RNA-seq Analysis Pipeline
# ============================================================

# STEP 1: Load packages
# --------------------------------------------------------
# These are like "apps" that give R extra abilities

library(DESeq2)       # For differential expression
library(ggplot2)      # For plotting
library(pheatmap)     # For heatmaps
library(dplyr)        # For data manipulation
library(tibble)       # For data frames
library(stringr)      # For text manipulation
library(dplyr)
select <- dplyr::select  # Force dplyr's select to be default

message("✓ Packages loaded successfully!")

# STEP 2: Set working directory
# --------------------------------------------------------
# Point this to your local syn2580853/RNA_Seq folder after downloading
# count files from Synapse into Collapsed_Counts/.
setwd("Your_Working_Directory/RNA_Seq")

message("Working directory: ", getwd())

# STEP 3: Create output folders
# --------------------------------------------------------
OUT_DIR <- "RESULTS_syn2580853"
GSEA_DIR <- "GSEA"
EXTRAS_DIR <- "extras"
dir.create(file.path(OUT_DIR, "figs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(GSEA_DIR, "figs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(GSEA_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(EXTRAS_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
message("✓ Output folders created")

# STEP 4: Find your counts files
# --------------------------------------------------------
count_files <- list.files("Collapsed_Counts", pattern = "_counts\\.txt$", full.names = TRUE)
message("Found ", length(count_files), " count files")

# Check: Did it find 18 files?
if (length(count_files) == 0) {
  stop("ERROR: No count files found! Check that Collapsed_Counts/ contains the _counts.txt files")
}

# Show the files found
print(basename(count_files))

# STEP 5: Read and merge all counts into one table
# --------------------------------------------------------
message("Reading count files...")

# Function to read one counts file
read_counts <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t", 
                   col.names = c("gene", "count"),
                   stringsAsFactors = FALSE)
  return(df)
}

# Get sample names from filenames (e.g., "DE01D01" from "Trim_DE01D01_counts.txt")
sample_names <- basename(count_files) %>%
  str_remove("^Trim_") %>%
  str_remove("_counts\\.txt$")

message("Sample names: ", paste(sample_names, collapse = ", "))

# Read all files into a list
counts_list <- lapply(count_files, read_counts)
names(counts_list) <- sample_names

# Merge all into one data frame
counts_merged <- counts_list[[1]] %>% rename(!!sample_names[1] := count)

for (i in 2:length(counts_list)) {
  temp <- counts_list[[i]] %>% rename(!!sample_names[i] := count)
  counts_merged <- full_join(counts_merged, temp, by = "gene")
}

# Convert to matrix format
rownames(counts_merged) <- counts_merged$gene
counts_mat <- as.matrix(counts_merged[, -1])
mode(counts_mat) <- "integer"
counts_mat[is.na(counts_mat)] <- 0

message("✓ Count matrix created: ", nrow(counts_mat), " genes x ", ncol(counts_mat), " samples")

# STEP 6: Create metadata from sample names
# --------------------------------------------------------
# Sample naming: DE01D01
#   DE = ELAV Control (E), DT = Tau (T)
#   01/10/20 = Day 1, 10, 20
#   D01/D02/D03 = replicate number

meta <- data.frame(
  sample = sample_names,
  genotype = ifelse(str_detect(sample_names, "^DE"), "Control", "Tau_WT"),
  day = case_when(
    str_detect(sample_names, "^..01") ~ "Day01",
    str_detect(sample_names, "^..10") ~ "Day10", 
    str_detect(sample_names, "^..20") ~ "Day20"
  ),
  replicate = str_extract(sample_names, "D0[0-9]$"),
  row.names = sample_names
)

# Set factor levels (Control is reference)
meta$genotype <- factor(meta$genotype, levels = c("Control", "Tau_WT"))
meta$day <- factor(meta$day, levels = c("Day01", "Day10", "Day20"))

# Display metadata
message("✓ Metadata created:")
print(meta)

# Save metadata
write.csv(meta, file.path(OUT_DIR, "tables", "sample_metadata.csv"), row.names = FALSE)
write.csv(meta, file.path("Collapsed_Counts", "metadata_from_filenames.csv"), row.names = FALSE)

# STEP 7: Quality Control - Library Sizes
# --------------------------------------------------------
message("Running QC...")

lib_sizes <- colSums(counts_mat)

# Create bar plot
png(file.path(OUT_DIR, "figs", "01_library_sizes.png"), width = 1000, height = 600, res = 150)
par(mar = c(8, 4, 3, 1))
barplot(lib_sizes / 1e6, 
        las = 2, 
        ylab = "Library Size (millions)",
        main = "Library Sizes per Sample",
        col = ifelse(meta$genotype == "Control", "steelblue", "salmon"),
        cex.names = 0.8)
legend("topright", legend = c("Control", "Tau_WT"), fill = c("steelblue", "salmon"))
dev.off()

message("✓ Library size plot saved")

# STEP 8: Filter low-count genes
# --------------------------------------------------------
# Keep genes with at least 10 counts in at least 20% of samples
min_samples <- ceiling(ncol(counts_mat) * 0.2)
keep <- rowSums(counts_mat >= 10) >= min_samples
counts_filt <- counts_mat[keep, ]

message("Genes before filtering: ", nrow(counts_mat))
message("Genes after filtering: ", nrow(counts_filt))
message("Removed ", nrow(counts_mat) - nrow(counts_filt), " low-count genes")

# STEP 9: Create DESeq2 object
# --------------------------------------------------------
message("Creating DESeq2 dataset...")

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = meta,
  design = ~ day + genotype  # Test genotype effect, accounting for day
)

# STEP 10: Run differential expression analysis
# --------------------------------------------------------
message("Running DESeq2 (this may take a minute)...")

dds <- DESeq(dds)

message("✓ DESeq2 complete!")

# STEP 11: Transform data for visualization
# --------------------------------------------------------
vsd <- vst(dds, blind = TRUE)

# STEP 12: PCA Plot
# --------------------------------------------------------
pca_data <- plotPCA(vsd, intgroup = c("genotype", "day"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

png(file.path(OUT_DIR, "figs", "02_PCA_plot.png"), width = 1000, height = 800, res = 150)
ggplot(pca_data, aes(x = PC1, y = PC2, color = genotype, shape = day)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("Control" = "steelblue", "Tau_WT" = "salmon")) +
  theme_minimal(base_size = 14) +
  ggtitle("PCA: Samples by Genotype and Day") +
  theme(legend.position = "right")
dev.off()

message("✓ PCA plot saved")

# STEP 13: Get differential expression results
# --------------------------------------------------------
message("Extracting results: Tau vs Control...")

res <- results(dds, 
               contrast = c("genotype", "Tau_WT", "Control"),
               alpha = 0.05)

# View summary
summary(res)

# Convert to data frame and sort by adjusted p-value
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# Save full results
write.csv(res_df, 
          file.path(OUT_DIR, "tables", "DE_Tau_vs_Control.csv"),
          row.names = FALSE)

# Count significant genes
n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
n_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)

message("✓ Significant DE genes (padj < 0.05): ", n_sig)
message("  - Upregulated in Tau: ", n_up)
message("  - Downregulated in Tau: ", n_down)

# STEP 14: Volcano Plot
# --------------------------------------------------------
res_plot <- res_df %>%
  mutate(
    significance = case_when(
      is.na(padj) ~ "NS",
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

png(file.path(OUT_DIR, "figs", "03_volcano_plot.png"), width = 1000, height = 800, res = 150)
ggplot(res_plot, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray70")) +
  theme_minimal(base_size = 14) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  ggtitle("Volcano Plot: Tau vs Control") +
  xlab("log2 Fold Change") +
  ylab("-log10(p-value)")
dev.off()

message("✓ Volcano plot saved")

# STEP 15: Top DE genes heatmap
# --------------------------------------------------------
top_genes <- res_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

if (length(top_genes) >= 5) {
  mat <- assay(vsd)[top_genes, ]
  mat_z <- t(scale(t(mat)))  # Z-score normalize
  
  anno_col <- data.frame(
    Genotype = meta$genotype,
    Day = meta$day,
    row.names = rownames(meta)
  )
  
  png(file.path(OUT_DIR, "figs", "04_heatmap_top50.png"), width = 1000, height = 1200, res = 150)
  pheatmap(mat_z,
           annotation_col = anno_col,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           show_rownames = TRUE,
           fontsize_row = 8,
           main = "Top 50 DE Genes")
  dev.off()
  
  message("✓ Heatmap saved")
} else {
  message("Not enough significant genes for heatmap")
}

# ============================================================
# DONE!
# ============================================================
message("\n========================================")
message("ANALYSIS COMPLETE!")
message("========================================")
message("Results saved to: ", OUT_DIR)
message("\nOutput files:")
message("  - tables/sample_metadata.csv")
message("  - tables/DE_Tau_vs_Control.csv")
message("  - figs/01_library_sizes.png")
message("  - figs/02_PCA_plot.png")
message("  - figs/03_volcano_plot.png")
message("  - figs/04_heatmap_top50.png")
message("\nTop 10 DE genes:")
print(head(res_df, 10))

# ============================================================
# PART 2: TIMEPOINT-SPECIFIC ANALYSIS & ENHANCED VISUALIZATION
# Following OSD-514 Procedure
# ============================================================

# Install additional packages if needed
if (!requireNamespace("apeglm", quietly = TRUE)) {
  BiocManager::install("apeglm", ask = FALSE)
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
if (!requireNamespace("enrichR", quietly = TRUE)) {
  install.packages("enrichR")
}

library(apeglm)
library(ggrepel)

# --------------------------------------------------------
# STEP 16: LFC Shrinkage (apeglm) for better fold change estimates
# --------------------------------------------------------
message("Applying LFC shrinkage...")

# Get the coefficient name for Tau vs Control
resultsNames(dds)

res_shrunk <- lfcShrink(dds, 
                        coef = "genotype_Tau_WT_vs_Control", 
                        type = "apeglm")

# Save shrunk results
res_shrunk_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

write.csv(res_shrunk_df, 
          file.path(OUT_DIR, "tables", "DE_Tau_vs_Control_shrunk.csv"),
          row.names = FALSE)

message("✓ LFC shrinkage complete")

# --------------------------------------------------------
# STEP 17: Timepoint-Specific Comparisons
# --------------------------------------------------------
message("Running timepoint-specific analyses...")

# Create a combined factor for genotype:day interaction
meta$group <- factor(paste(meta$genotype, meta$day, sep = "_"))

# Rebuild DESeq2 with interaction design
dds_time <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = meta,
  design = ~ group
)
dds_time <- DESeq(dds_time)

# Function to run contrast and save results
run_contrast <- function(dds, contrast_name, num, denom) {
  res <- results(dds, contrast = c("group", num, denom), alpha = 0.05)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)
  
  # Save full results
  write.csv(res_df, 
            file.path(OUT_DIR, "tables", paste0("DE_", contrast_name, ".csv")),
            row.names = FALSE)
  
  # Count significant
  n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  n_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)
  n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)
  
  message(sprintf("%s: %d significant (↑%d, ↓%d)", 
                  contrast_name, n_sig, n_up, n_down))
  
  return(list(res = res, res_df = res_df))
}

# Run timepoint-specific contrasts
res_Day01 <- run_contrast(dds_time, "Tau_vs_Control_Day01", 
                          "Tau_WT_Day01", "Control_Day01")
res_Day10 <- run_contrast(dds_time, "Tau_vs_Control_Day10", 
                          "Tau_WT_Day10", "Control_Day10")
res_Day20 <- run_contrast(dds_time, "Tau_vs_Control_Day20", 
                          "Tau_WT_Day20", "Control_Day20")

message("✓ Timepoint-specific analyses complete")

# --------------------------------------------------------
# STEP 18: Enhanced Volcano Plots (with gene labels)
# --------------------------------------------------------
message("Creating enhanced volcano plots...")

save_volcano_enhanced <- function(res_df, tag, top_n_labels = 15, 
                                  thr_p = 0.05, thr_fc = 1) {
  
  vdf <- res_df %>%
    mutate(
      group = case_when(
        is.na(padj) ~ "NS",
        padj < thr_p & log2FoldChange >= thr_fc ~ "Up (sig)",
        padj < thr_p & log2FoldChange <= -thr_fc ~ "Down (sig)",
        TRUE ~ "NS"
      ),
      group = factor(group, levels = c("Down (sig)", "NS", "Up (sig)"))
    )
  
  # Select top genes to label
  top_genes <- vdf %>%
    filter(padj < thr_p) %>%
    arrange(padj) %>%
    head(top_n_labels) %>%
    pull(gene)
  
  p <- ggplot(vdf, aes(x = log2FoldChange, y = -log10(pvalue), color = group)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-thr_fc, thr_fc), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(thr_p), linetype = "dashed", color = "gray40") +
    geom_text_repel(
      data = subset(vdf, gene %in% top_genes),
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5
    ) +
    scale_color_manual(values = c("Down (sig)" = "#2C7BB6", 
                                  "NS" = "grey70", 
                                  "Up (sig)" = "#D7191C")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top") +
    labs(
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      title = paste("Volcano:", gsub("_", " ", tag)),
      subtitle = "Dashed lines: |LFC| = 1, p = 0.05"
    )
  
  ggsave(file.path(OUT_DIR, "figs", paste0("volcano_", tag, ".png")),
         p, width = 10, height = 8, dpi = 200)
  
  message("  Saved: volcano_", tag, ".png")
}

# Create volcano plots for each comparison
save_volcano_enhanced(res_df, "Tau_vs_Control_Overall")
save_volcano_enhanced(res_Day01$res_df, "Tau_vs_Control_Day01")
save_volcano_enhanced(res_Day10$res_df, "Tau_vs_Control_Day10")
save_volcano_enhanced(res_Day20$res_df, "Tau_vs_Control_Day20")

message("✓ Volcano plots saved")

# --------------------------------------------------------
# STEP 19: Enhanced Heatmaps (with sample annotations)
# --------------------------------------------------------
message("Creating enhanced heatmaps...")

save_heatmap_enhanced <- function(res_df, tag, topN = 40, cap_z = 2.5) {
  
  # Get top genes
  top_genes <- res_df %>%
    filter(padj < 0.05) %>%
    arrange(padj, -abs(log2FoldChange)) %>%
    head(topN) %>%
    pull(gene)
  
  if (length(top_genes) < 5) {
    message("  Not enough significant genes for heatmap: ", tag)
    return(invisible(NULL))
  }
  
  # Get expression matrix
  mat <- assay(vsd)[top_genes, ]
  
  # Z-score normalize and cap
  mat_z <- t(scale(t(mat)))
  mat_z[mat_z > cap_z] <- cap_z
  mat_z[mat_z < -cap_z] <- -cap_z
  
  # Annotation
  ann <- data.frame(
    Genotype = meta$genotype,
    Day = meta$day,
    row.names = rownames(meta)
  )
  
  # Order columns by genotype then day
  col_order <- order(meta$genotype, meta$day)
  mat_z <- mat_z[, col_order]
  ann <- ann[col_order, , drop = FALSE]
  
  # Create short labels for columns
  short_labels <- paste0(
    ifelse(meta$genotype == "Control", "C", "T"),
    "_",
    gsub("Day", "D", meta$day),
    "_",
    meta$replicate
  )
  names(short_labels) <- rownames(meta)
  
  # Color palette
  ann_colors <- list(
    Genotype = c(Control = "steelblue", Tau_WT = "salmon"),
    Day = c(Day01 = "#fee8c8", Day10 = "#fdbb84", Day20 = "#e34a33")
  )
  
  png(file.path(OUT_DIR, "figs", paste0("heatmap_", tag, ".png")),
      width = 1400, height = 1200, res = 150)
  
  pheatmap(mat_z,
           annotation_col = ann,
           annotation_colors = ann_colors,
           labels_col = short_labels[colnames(mat_z)],
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_rownames = TRUE,
           fontsize_row = 8,
           fontsize_col = 9,
           border_color = NA,
           main = paste0("Top ", length(top_genes), " DE Genes: ", 
                         gsub("_", " ", tag), " (Z-score, cap ±", cap_z, ")"))
  
  dev.off()
  
  message("  Saved: heatmap_", tag, ".png")
}

# Create heatmaps
save_heatmap_enhanced(res_df, "Tau_vs_Control_Overall", topN = 50)
save_heatmap_enhanced(res_Day01$res_df, "Tau_vs_Control_Day01", topN = 40)
save_heatmap_enhanced(res_Day10$res_df, "Tau_vs_Control_Day10", topN = 40)
save_heatmap_enhanced(res_Day20$res_df, "Tau_vs_Control_Day20", topN = 40)

message("✓ Heatmaps saved")

# --------------------------------------------------------
# STEP 20: Mitochondrial Gene Analysis
# --------------------------------------------------------
message("Running mitochondrial gene analysis...")

# Install org.Dm.eg.db if needed
if (!requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Dm.eg.db", ask = FALSE)
}
library(org.Dm.eg.db)
library(AnnotationDbi)
library(GO.db)

# Curated mitochondrial GO terms (from OSD-514 procedure)
mito_go_bp <- c(
  "GO:0006119",  # oxidative phosphorylation
  "GO:0022900",  # electron transport chain
  "GO:0006099",  # tricarboxylic acid (TCA) cycle
  "GO:0006635",  # fatty acid beta-oxidation
  "GO:0000422"   # mitophagy
)

mito_go_cc <- c(
  "GO:0005739",  # mitochondrion
  "GO:0005743",  # mitochondrial inner membrane
  "GO:0005747",  # respiratory chain complex I
  "GO:0005753",  # ATP synthase complex
  "GO:0005759"   # mitochondrial matrix
)

# Get all genes annotated with mito GO terms
# Note: Your data uses gene symbols, so we'll query by symbol
all_mito_go <- c(mito_go_bp, mito_go_cc)

# Get genes for each GO term
mito_genes_list <- lapply(all_mito_go, function(go_id) {
  tryCatch({
    genes <- AnnotationDbi::select(org.Dm.eg.db, 
                                   keys = go_id, 
                                   keytype = "GO",
                                   columns = "SYMBOL")
    unique(na.omit(genes$SYMBOL))
  }, error = function(e) character(0))
})

# Combine all mito genes
mito_genes_curated <- unique(unlist(mito_genes_list))
message("Curated mitochondrial gene set: ", length(mito_genes_curated), " genes")

# Find overlap with our DE genes
mito_in_data <- intersect(mito_genes_curated, rownames(counts_filt))
message("Mito genes in our dataset: ", length(mito_in_data))

# Check which are significant
mito_sig <- res_df %>%
  filter(gene %in% mito_in_data & padj < 0.05)

message("Significantly DE mito genes: ", nrow(mito_sig))

# Save mito gene results
mito_results <- res_df %>%
  filter(gene %in% mito_in_data) %>%
  arrange(padj)

write.csv(mito_results,
          file.path(OUT_DIR, "tables", "DE_mitochondrial_genes.csv"),
          row.names = FALSE)

# --------------------------------------------------------
# STEP 21: Mitochondrial Gene Volcano Plot
# --------------------------------------------------------
save_volcano_mito <- function(res_df, mito_genes, tag, thr_p = 0.05, thr_fc = 1) {
  
  vdf <- res_df %>%
    mutate(
      is_mito = gene %in% mito_genes,
      group = case_when(
        is.na(padj) ~ "NS",
        padj < thr_p & log2FoldChange >= thr_fc ~ "Up (sig)",
        padj < thr_p & log2FoldChange <= -thr_fc ~ "Down (sig)",
        TRUE ~ "NS"
      ),
      group = factor(group, levels = c("Down (sig)", "NS", "Up (sig)"))
    )
  
  # Get mito genes to label (significant ones)
  mito_to_label <- vdf %>%
    filter(is_mito & padj < 0.1) %>%
    arrange(padj) %>%
    head(20) %>%
    pull(gene)
  
  p <- ggplot(vdf, aes(x = log2FoldChange, y = -log10(pvalue), color = group)) +
    geom_point(alpha = 0.5, size = 1.5) +
    # Highlight mito genes with black outline
    geom_point(data = subset(vdf, is_mito),
               shape = 21, stroke = 0.8, size = 2.5, 
               aes(fill = group), color = "black") +
    geom_vline(xintercept = c(-thr_fc, thr_fc), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(thr_p), linetype = "dashed", color = "gray40") +
    geom_text_repel(
      data = subset(vdf, gene %in% mito_to_label),
      aes(label = gene),
      size = 3,
      max.overlaps = 30,
      box.padding = 0.4
    ) +
    scale_color_manual(values = c("Down (sig)" = "#2C7BB6", 
                                  "NS" = "grey70", 
                                  "Up (sig)" = "#D7191C")) +
    scale_fill_manual(values = c("Down (sig)" = "#2C7BB6", 
                                 "NS" = "grey70", 
                                 "Up (sig)" = "#D7191C")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top") +
    labs(
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      title = paste("Volcano (Mito Focus):", gsub("_", " ", tag)),
      subtitle = "Black-outlined points = mitochondrial genes"
    )
  
  ggsave(file.path(OUT_DIR, "figs", paste0("volcano_mito_", tag, ".png")),
         p, width = 10, height = 8, dpi = 200)
  
  message("  Saved: volcano_mito_", tag, ".png")
}

save_volcano_mito(res_df, mito_in_data, "Tau_vs_Control")

# --------------------------------------------------------
# STEP 22: Mitochondrial Gene Heatmap
# --------------------------------------------------------
if (length(mito_in_data) >= 10) {
  # Get DE mito genes (relax threshold slightly)
  mito_de <- res_df %>%
    filter(gene %in% mito_in_data & padj < 0.1) %>%
    arrange(padj) %>%
    head(50) %>%
    pull(gene)
  
  if (length(mito_de) >= 5) {
    mat <- assay(vsd)[mito_de, ]
    mat_z <- t(scale(t(mat)))
    mat_z[mat_z > 2.5] <- 2.5
    mat_z[mat_z < -2.5] <- -2.5
    
    ann <- data.frame(
      Genotype = meta$genotype,
      Day = meta$day,
      row.names = rownames(meta)
    )
    
    col_order <- order(meta$genotype, meta$day)
    mat_z <- mat_z[, col_order]
    ann <- ann[col_order, , drop = FALSE]
    
    ann_colors <- list(
      Genotype = c(Control = "steelblue", Tau_WT = "salmon"),
      Day = c(Day01 = "#fee8c8", Day10 = "#fdbb84", Day20 = "#e34a33")
    )
    
    png(file.path(OUT_DIR, "figs", "heatmap_mitochondrial_genes.png"),
        width = 1400, height = 1200, res = 150)
    
    pheatmap(mat_z,
             annotation_col = ann,
             annotation_colors = ann_colors,
             cluster_cols = FALSE,
             cluster_rows = TRUE,
             show_rownames = TRUE,
             fontsize_row = 8,
             border_color = NA,
             main = paste0("Mitochondrial DE Genes (", length(mito_de), " genes)"))
    
    dev.off()
    
    message("✓ Mitochondrial heatmap saved")
  }
}

message("\n✓ Mitochondrial analysis complete")

# --------------------------------------------------------
# STEP 23: GO Enrichment Analysis (using enrichR)
# --------------------------------------------------------
message("Running GO enrichment analysis...")

library(enrichR)

# Set up enrichR databases
dbs <- c("GO_Biological_Process_2021",
         "GO_Cellular_Component_2021", 
         "GO_Molecular_Function_2021",
         "KEGG_2021_Drosophila")

# Get significant gene lists
sig_up <- res_df %>% 
  filter(padj < 0.05 & log2FoldChange > 1) %>% 
  pull(gene)

sig_down <- res_df %>% 
  filter(padj < 0.05 & log2FoldChange < -1) %>% 
  pull(gene)

sig_all <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(gene)

message("Upregulated genes: ", length(sig_up))
message("Downregulated genes: ", length(sig_down))

# Run enrichment
if (length(sig_all) >= 10) {
  
  enrich_all <- enrichr(sig_all, dbs)
  
  # Save results
  for (db_name in names(enrich_all)) {
    write.csv(enrich_all[[db_name]],
              file.path(EXTRAS_DIR, "tables",
                        paste0("GO_", gsub(" ", "_", db_name), "_all_DE.csv")),
              row.names = FALSE)
  }
  
  # Plot top GO Biological Process terms
  if (nrow(enrich_all[["GO_Biological_Process_2021"]]) > 0) {
    top_bp <- enrich_all[["GO_Biological_Process_2021"]] %>%
      filter(Adjusted.P.value < 0.05) %>%
      arrange(Adjusted.P.value) %>%
      head(20) %>%
      mutate(Term = substr(Term, 1, 50))  # Truncate long names
    
    if (nrow(top_bp) > 0) {
      p <- ggplot(top_bp, aes(x = -log10(Adjusted.P.value), 
                              y = reorder(Term, -log10(Adjusted.P.value)))) +
        geom_bar(stat = "identity", fill = "steelblue") +
        theme_minimal(base_size = 11) +
        labs(x = "-log10(FDR)", 
             y = "", 
             title = "GO Biological Process Enrichment",
             subtitle = "Tau vs Control DE Genes") +
        theme(axis.text.y = element_text(size = 9))
      
      ggsave(file.path(OUT_DIR, "figs", "GO_BP_enrichment.png"),
             p, width = 10, height = 8, dpi = 200)
      
      message("✓ GO enrichment plot saved")
    }
  }
  
  # Identify mito-related enriched terms
  mito_terms <- enrich_all[["GO_Biological_Process_2021"]] %>%
    filter(grepl("mitochond|respiratory|oxidative|electron|ATP", 
                 Term, ignore.case = TRUE))
  
  if (nrow(mito_terms) > 0) {
    write.csv(mito_terms,
              file.path(EXTRAS_DIR, "tables", "GO_mitochondrial_terms.csv"),
              row.names = FALSE)
    message("Found ", nrow(mito_terms), " mitochondria-related GO terms")
  }
  
} else {
  message("Not enough DE genes for enrichment analysis")
}

# ============================================================
# ANALYSIS COMPLETE
# ============================================================
message("\n========================================")
message("FULL ANALYSIS COMPLETE!")
message("========================================")
message("\nOutput files in: ", OUT_DIR)
message("\nTables generated:")
message("  - DE_Tau_vs_Control.csv (overall)")
message("  - DE_Tau_vs_Control_shrunk.csv (with LFC shrinkage)")
message("  - DE_Tau_vs_Control_Day01/10/20.csv (timepoint-specific)")
message("  - DE_mitochondrial_genes.csv")
message("  - GO enrichment tables")
message("\nFigures generated:")
message("  - Volcano plots (overall + timepoint-specific + mito focus)")
message("  - Heatmaps (overall + timepoint-specific + mito)")
message("  - GO enrichment barplot")

# Print summary
message("\n===== SUMMARY =====")
message("Total DE genes (padj < 0.05): ", sum(res_df$padj < 0.05, na.rm = TRUE))
message("  Day01: ", sum(res_Day01$res_df$padj < 0.05, na.rm = TRUE))
message("  Day10: ", sum(res_Day10$res_df$padj < 0.05, na.rm = TRUE))
message("  Day20: ", sum(res_Day20$res_df$padj < 0.05, na.rm = TRUE))
message("Mito genes in dataset: ", length(mito_in_data))
message("Significant mito genes: ", nrow(mito_sig))

# ============================================================
# PART 3: GSEA (Gene Set Enrichment Analysis) with fgsea
# Following OSD-514 Procedure
# ============================================================

message("Starting GSEA analysis...")

# Install fgsea if needed
if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", ask = FALSE)
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

library(fgsea)
library(data.table)

# --------------------------------------------------------
# STEP 24: Build Gene Ranks
# --------------------------------------------------------
# GSEA uses ALL genes ranked by a metric (not just significant ones)
# We'll use: sign(log2FC) * -log10(pvalue) as the ranking metric

build_ranks <- function(res_df) {
  # Remove NAs
  df <- res_df %>%
    filter(!is.na(pvalue) & !is.na(log2FoldChange) & pvalue > 0)
  
  # Create ranking metric: signed -log10(p-value)
  ranks <- sign(df$log2FoldChange) * -log10(df$pvalue)
  names(ranks) <- df$gene
  
  # Sort by rank
  ranks <- sort(ranks, decreasing = TRUE)
  
  return(ranks)
}

# Build ranks for overall comparison
ranks_overall <- build_ranks(res_df)
message("Ranked ", length(ranks_overall), " genes for GSEA")

# --------------------------------------------------------
# STEP 25: Build GO Biological Process Gene Sets
# --------------------------------------------------------
message("Building GO BP gene sets...")

# Get all GO BP terms for genes in our dataset
gene_list <- names(ranks_overall)

# Query GO annotations for our genes
go_annotations <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = gene_list,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GO", "ONTOLOGY")
)

# Filter to Biological Process only - WITH dplyr:: prefix
go_bp <- go_annotations %>%
  dplyr::filter(ONTOLOGY == "BP") %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::select(GO, SYMBOL) %>%
  dplyr::distinct()

# Create pathway list (GO term -> gene vector)
pathways_bp <- split(go_bp$SYMBOL, go_bp$GO)

# Filter pathways by size (10-500 genes)
pathway_sizes <- sapply(pathways_bp, length)
pathways_bp <- pathways_bp[pathway_sizes >= 10 & pathway_sizes <= 500]

message("Built ", length(pathways_bp), " GO BP gene sets (size 10-500)")

# Get GO term names for labeling
go_term_names <- AnnotationDbi::select(
  GO.db,
  keys = names(pathways_bp),
  keytype = "GOID",
  columns = "TERM"
)
term_map <- setNames(go_term_names$TERM, go_term_names$GOID)

# --------------------------------------------------------
# STEP 26: Run fgsea
# --------------------------------------------------------
message("Running fgsea (this may take a minute)...")

set.seed(42)
fgsea_results <- fgsea(
  pathways = pathways_bp,
  stats = ranks_overall,
  minSize = 10,
  maxSize = 500,
  nPermSimple = 10000
)

# Add term names
fgsea_results$term <- term_map[fgsea_results$pathway]

# Sort by adjusted p-value
fgsea_results <- fgsea_results[order(padj)]

message("GSEA complete: ", sum(fgsea_results$padj < 0.05, na.rm = TRUE), 
        " significant pathways (FDR < 0.05)")

# --------------------------------------------------------
# STEP 27: Save GSEA Results
# --------------------------------------------------------

# Save full results (collapse leadingEdge to string)
fgsea_export <- fgsea_results %>%
  as.data.frame() %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write.csv(fgsea_export,
          file.path(GSEA_DIR, "tables", "fgsea_GO_BP_Tau_vs_Control.csv"),
          row.names = FALSE)

# Save significant results
fgsea_sig <- fgsea_export %>% filter(padj < 0.05)
write.csv(fgsea_sig,
          file.path(GSEA_DIR, "tables", "fgsea_GO_BP_significant.csv"),
          row.names = FALSE)

message("Saved GSEA tables to: ", file.path(GSEA_DIR, "tables"))

# --------------------------------------------------------
# STEP 28: GSEA Barplot (Top Pathways by NES)
# --------------------------------------------------------
message("Creating GSEA plots...")

# Get top positive and negative NES pathways
top_up <- fgsea_results %>%
  filter(padj < 0.1 & NES > 0) %>%
  arrange(desc(NES)) %>%
  head(10)

top_down <- fgsea_results %>%
  filter(padj < 0.1 & NES < 0) %>%
  arrange(NES) %>%
  head(10)

top_pathways <- rbind(top_up, top_down) %>%
  mutate(term_short = ifelse(nchar(term) > 50, 
                             paste0(substr(term, 1, 47), "..."), 
                             term))

if (nrow(top_pathways) > 0) {
  p <- ggplot(top_pathways, aes(x = reorder(term_short, NES), y = NES, 
                                fill = -log10(padj))) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal(base_size = 11) +
    labs(
      x = "",
      y = "Normalized Enrichment Score (NES)",
      fill = "-log10(FDR)",
      title = "GSEA: GO Biological Process",
      subtitle = "Tau vs Control - Top Enriched Pathways"
    ) +
    theme(axis.text.y = element_text(size = 9))
  
  ggsave(file.path(GSEA_DIR, "figs", "fgsea_barplot_top_pathways.png"),
         p, width = 12, height = 8, dpi = 200)
  
  message("✓ Saved: fgsea_barplot_top_pathways.png")
}

# --------------------------------------------------------
# STEP 29: Identify Mito-Related GSEA Pathways
# --------------------------------------------------------
message("Identifying mitochondrial pathways in GSEA results...")

# Pattern to match mito-related terms
mito_pattern <- paste(
  "mitochond", "oxidative phosph", "electron transport",
  "respiratory chain", "ATP synthase", "tricarboxylic", "TCA",
  "beta-oxid", "mitophagy", "cristae",
  sep = "|"
)

fgsea_mito <- fgsea_results %>%
  filter(grepl(mito_pattern, term, ignore.case = TRUE))

message("Found ", nrow(fgsea_mito), " mitochondria-related pathways")
message("  Significant (FDR < 0.05): ", sum(fgsea_mito$padj < 0.05, na.rm = TRUE))

# Save mito pathways
if (nrow(fgsea_mito) > 0) {
  fgsea_mito_export <- fgsea_mito %>%
    as.data.frame() %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
  
  write.csv(fgsea_mito_export,
            file.path(GSEA_DIR, "tables", "fgsea_mitochondrial_pathways.csv"),
            row.names = FALSE)
  
  # Print top mito pathways
  message("\nTop Mitochondrial Pathways:")
  print(fgsea_mito %>% 
          as.data.frame() %>%
          dplyr::select(term, NES, padj) %>% 
          head(10))
}

# --------------------------------------------------------
# STEP 30: Enrichment Curves for Top Mito Pathways
# --------------------------------------------------------

# Plot enrichment curves for top mito-related pathways
if (nrow(fgsea_mito) > 0) {
  
  # Get top 4 mito pathways (by significance)
  top_mito <- fgsea_mito %>%
    arrange(padj) %>%
    head(4) %>%
    pull(pathway)
  
  for (pw in top_mito) {
    pw_name <- term_map[pw]
    if (is.na(pw_name)) pw_name <- pw
    
    # Create enrichment plot
    p <- plotEnrichment(pathways_bp[[pw]], ranks_overall) +
      ggtitle(paste0(pw_name, "\nNES = ", 
                     round(fgsea_results[pathway == pw, NES], 2),
                     ", FDR = ", 
                     signif(fgsea_results[pathway == pw, padj], 2)))
    
    # Clean filename
    fname <- gsub("[^A-Za-z0-9]", "_", substr(pw_name, 1, 40))
    
    ggsave(file.path(GSEA_DIR, "figs", paste0("enrichment_curve_", fname, ".png")),
           p, width = 8, height = 5, dpi = 200)
  }
  
  message("✓ Saved enrichment curves for top mito pathways")
}

# --------------------------------------------------------
# STEP 31: GSEA Table Plot (fgsea style)
# --------------------------------------------------------

# Create table plot for top pathways
if (nrow(top_pathways) > 0) {
  top_pw_ids <- top_pathways$pathway
  
  png(file.path(GSEA_DIR, "figs", "fgsea_table_top_pathways.png"),
      width = 1600, height = 900, res = 150)
  
  plotGseaTable(
    pathways_bp[top_pw_ids],
    ranks_overall,
    fgsea_results[pathway %in% top_pw_ids],
    gseaParam = 0.5
  )
  
  dev.off()
  
  message("✓ Saved: fgsea_table_top_pathways.png")
}

# --------------------------------------------------------
# STEP 32: Run GSEA for Each Timepoint
# --------------------------------------------------------
message("\nRunning GSEA for each timepoint...")

run_gsea_timepoint <- function(res_df, tag) {
  ranks <- build_ranks(res_df)
  
  if (length(ranks) < 100) {
    message("  ", tag, ": Too few ranked genes, skipping")
    return(NULL)
  }
  
  set.seed(42)
  fg <- fgsea(
    pathways = pathways_bp,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    nPermSimple = 10000
  )
  
  fg$term <- term_map[fg$pathway]
  fg <- fg[order(padj)]
  
  n_sig <- sum(fg$padj < 0.05, na.rm = TRUE)
  message("  ", tag, ": ", n_sig, " significant pathways")
  
  # Save results
  fg_export <- fg %>%
    as.data.frame() %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
  
  write.csv(fg_export,
            file.path(GSEA_DIR, "tables", paste0("fgsea_GO_BP_", tag, ".csv")),
            row.names = FALSE)
  
  return(fg)
}

fgsea_Day01 <- run_gsea_timepoint(res_Day01$res_df, "Day01")
fgsea_Day10 <- run_gsea_timepoint(res_Day10$res_df, "Day10")
fgsea_Day20 <- run_gsea_timepoint(res_Day20$res_df, "Day20")

# ============================================================
# GSEA COMPLETE
# ============================================================
message("\n========================================")
message("GSEA ANALYSIS COMPLETE!")
message("========================================")
message("\nOutput in: ", GSEA_DIR)
message("\nTables:")
message("  - fgsea_GO_BP_Tau_vs_Control.csv (all pathways)")
message("  - fgsea_GO_BP_significant.csv (FDR < 0.05)")
message("  - fgsea_mitochondrial_pathways.csv")
message("  - fgsea_GO_BP_Day01/10/20.csv (timepoint-specific)")
message("\nFigures:")
message("  - fgsea_barplot_top_pathways.png")
message("  - fgsea_table_top_pathways.png")
message("  - enrichment_curve_*.png (mito pathways)")

message("\n===== GSEA SUMMARY =====")
message("Overall Tau vs Control:")
message("  Significant pathways (FDR < 0.05): ", 
        sum(fgsea_results$padj < 0.05, na.rm = TRUE))
message("  Mito-related pathways found: ", nrow(fgsea_mito))
message("  Mito-related significant: ", 
        sum(fgsea_mito$padj < 0.05, na.rm = TRUE))