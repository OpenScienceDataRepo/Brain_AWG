# ============================================================
# Emory Drosophila Tau Model - PROTEOMICS Analysis Pipeline
# Adapted from OSD-514 Proteomics Procedure
# ============================================================

# --------------------------------------------------------
# STEP 1: Load Packages
# --------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(limma)
library(pheatmap)
library(matrixStats)
library(tibble)

message("✓ Packages loaded")

# --------------------------------------------------------
# STEP 2: Set Up Directories
# --------------------------------------------------------
# Point this to your local syn2580853/Proteomics folder after downloading
# the Baylor TMT file and biospecimen metadata from Synapse.
setwd("Your_Working_Directory/Proteomics")

OUT_DIR <- "RESULTS_syn2580853"
EXTRAS_DIR <- "extras"
GSEA_DIR <- file.path(EXTRAS_DIR, "GSEA")
dir.create(file.path(OUT_DIR, "figs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(EXTRAS_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(GSEA_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(GSEA_DIR, "figs"), recursive = TRUE, showWarnings = FALSE)

message("✓ Output directory: ", OUT_DIR)

# --------------------------------------------------------
# STEP 3: Read Proteomics Data
# --------------------------------------------------------
message("Reading proteomics data...")

prot_raw <- read.table("TMT_all/baylor_proteomics_fly_proteinoutput.txt", 
                       header = TRUE, 
                       sep = "\t",
                       quote = "",
                       stringsAsFactors = FALSE)

message("Raw data: ", nrow(prot_raw), " proteins x ", ncol(prot_raw), " columns")

# --------------------------------------------------------
# STEP 4: Clean Protein IDs and Extract LFQ Intensities
# --------------------------------------------------------

lfq_cols <- grep("^LFQ", colnames(prot_raw), ignore.case = TRUE, value = TRUE)
message("Found ", length(lfq_cols), " LFQ intensity columns")

id_col <- if ("Protein.IDs" %in% colnames(prot_raw)) {
  "Protein.IDs"
} else {
  "Protein IDs"
}

# Extract expression matrix
expr_raw <- prot_raw[, lfq_cols]
rownames(expr_raw) <- prot_raw[[id_col]]

# Clean column names to sample numbers (01, 02, ...)
colnames(expr_raw) <- sub("^LFQ(?:[.]|\\s+)?intensity(?:[.]|\\s+)?", "", colnames(expr_raw), ignore.case = TRUE, perl = TRUE)
colnames(expr_raw) <- sprintf("%02d", as.integer(colnames(expr_raw)))

# Convert to numeric matrix
expr_raw <- as.matrix(expr_raw)
mode(expr_raw) <- "numeric"

# Replace 0 with NA (MaxQuant uses 0 for missing values)
expr_raw[expr_raw == 0] <- NA

message("Expression matrix: ", nrow(expr_raw), " proteins x ", ncol(expr_raw), " samples")

# --------------------------------------------------------
# STEP 5: Filter Contaminants and Reverse Hits
# --------------------------------------------------------
message("Filtering contaminants...")

# Remove contaminants (CON__) and reverse hits (REV__)
is_contaminant <- grepl("^CON__|^REV__", rownames(expr_raw))
expr_filt <- expr_raw[!is_contaminant, ]

message("Removed ", sum(is_contaminant), " contaminant/reverse proteins")
message("Remaining: ", nrow(expr_filt), " proteins")

# --------------------------------------------------------
# STEP 6: Build Metadata (CORRECTED)
# --------------------------------------------------------
message("Building metadata...")

# Read biospecimen metadata
bio_meta <- read.csv("EmoryDrosophilaTau_biospecimen_metadata.csv")

# Filter to proteomics samples (LC-MSMS)
prot_meta <- bio_meta %>%
  filter(assay == "LC-MSMS") %>%
  dplyr::select(specimenID, genotype, ageDays)

# Extract sample number from specimenID ("LFQ intensity 01" -> "01")
prot_meta$sample <- str_extract(prot_meta$specimenID, "\\d+$")

# Verify the mapping
message("Expression columns: ", paste(colnames(expr_filt), collapse = ", "))
message("Metadata samples:   ", paste(prot_meta$sample, collapse = ", "))

# Simplify genotype names
prot_meta <- prot_meta %>%
  mutate(
    genotype_simple = case_when(
      genotype == "Elav-GAL4/+" ~ "Control",
      genotype == "Elav-GAL4/+;UAS-TauWT/+" ~ "TauWT",
      genotype == "Elav-GAL4/+;UAS-TauR406W/+" ~ "TauR406W"
    ),
    day = paste0("Day", sprintf("%02d", ageDays)),
    group = paste(genotype_simple, day, sep = "_")
  )

# Set factor levels
prot_meta$genotype_simple <- factor(prot_meta$genotype_simple, 
                                    levels = c("Control", "TauWT", "TauR406W"))
prot_meta$day <- factor(prot_meta$day, 
                        levels = c("Day01", "Day10", "Day20"))
prot_meta$group <- factor(prot_meta$group)

# --- ALIGN: Keep only common samples ---
common_samples <- intersect(colnames(expr_filt), prot_meta$sample)
message("\nCommon samples found: ", length(common_samples))

if (length(common_samples) == 0) {
  stop("ERROR: No matching samples!")
}

# Subset and order
expr_filt <- expr_filt[, common_samples]
prot_meta <- prot_meta %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, colnames(expr_filt)))

rownames(prot_meta) <- prot_meta$sample

# --- VERIFY ---
stopifnot(all(colnames(expr_filt) == prot_meta$sample))
message("✓ Alignment verified! ", length(common_samples), " samples matched.")

message("✓ Metadata created")
print(table(prot_meta$genotype_simple, prot_meta$day))

# --------------------------------------------------------
# STEP 7: Log2 Transform and Filter
# --------------------------------------------------------
message("Log2 transforming and filtering...")

# Log2 transform
expr_log <- log2(expr_filt)

# Filter: keep proteins with values in at least 70% of samples
min_samples <- ceiling(0.7 * ncol(expr_log))
keep <- rowSums(!is.na(expr_log)) >= min_samples
expr_filt2 <- expr_log[keep, ]

# Also require some variance
expr_filt2 <- expr_filt2[rowVars(expr_filt2, na.rm = TRUE) > 0, ]

message("Proteins before filtering: ", nrow(expr_log))
message("Proteins after filtering: ", nrow(expr_filt2))

# --------------------------------------------------------
# STEP 8: Impute Missing Values (minimal imputation)
# --------------------------------------------------------
message("Imputing missing values...")

# For limma, impute missing values with row minimum - 1 (downshifted minimum)
expr_imputed <- expr_filt2
for (i in 1:nrow(expr_imputed)) {
  row_min <- min(expr_imputed[i, ], na.rm = TRUE)
  expr_imputed[i, is.na(expr_imputed[i, ])] <- row_min - 1
}

# --------------------------------------------------------
# STEP 9: Quality Control - PCA
# --------------------------------------------------------
message("Running PCA...")

pca <- prcomp(t(expr_imputed), center = TRUE, scale = FALSE)
pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(prot_meta, by = "sample")

percentVar <- round(100 * summary(pca)$importance[2, 1:2])

png(file.path(OUT_DIR, "figs", "01_PCA_proteomics.png"), 
    width = 1000, height = 800, res = 150)
ggplot(pca_df, aes(x = PC1, y = PC2, color = genotype_simple, shape = day)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Control" = "steelblue", 
                                "TauWT" = "salmon", 
                                "TauR406W" = "darkred")) +
  theme_minimal(base_size = 14) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    title = "PCA: Proteomics Samples",
    color = "Genotype",
    shape = "Day"
  )
dev.off()

message("✓ PCA plot saved")

# --------------------------------------------------------
# STEP 10: Differential Analysis with limma
# --------------------------------------------------------
message("Running limma differential analysis...")

# Create design matrix
design <- model.matrix(~ 0 + group, data = prot_meta)
colnames(design) <- gsub("group", "", colnames(design))

# Fit linear model
fit <- lmFit(expr_imputed, design)

# Define contrasts
contrast_matrix <- makeContrasts(
  # TauWT vs Control (overall across all days)
  TauWT_vs_Control_Day01 = TauWT_Day01 - Control_Day01,
  TauWT_vs_Control_Day10 = TauWT_Day10 - Control_Day10,
  TauWT_vs_Control_Day20 = TauWT_Day20 - Control_Day20,
  
  # TauR406W vs Control
  TauR406W_vs_Control_Day01 = TauR406W_Day01 - Control_Day01,
  TauR406W_vs_Control_Day10 = TauR406W_Day10 - Control_Day10,
  TauR406W_vs_Control_Day20 = TauR406W_Day20 - Control_Day20,
  
  # TauR406W vs TauWT (mutant vs wild-type tau)
  TauR406W_vs_TauWT_Day01 = TauR406W_Day01 - TauWT_Day01,
  TauR406W_vs_TauWT_Day10 = TauR406W_Day10 - TauWT_Day10,
  TauR406W_vs_TauWT_Day20 = TauR406W_Day20 - TauWT_Day20,
  
  levels = design
)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

message("✓ Limma analysis complete")

# --------------------------------------------------------
# STEP 11: Extract and Save Results
# --------------------------------------------------------
message("Extracting results...")

# Function to extract and save results
save_limma_results <- function(fit, coef_name) {
  tt <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("protein_id") %>%
    arrange(adj.P.Val)

  table_dir <- if (grepl("^TauR406W_vs_TauWT_", coef_name)) {
    file.path(EXTRAS_DIR, "tables")
  } else {
    file.path(OUT_DIR, "tables")
  }

  write.csv(tt,
            file.path(table_dir, paste0("DE_", coef_name, ".csv")),
            row.names = FALSE)
  
  n_sig <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0.5, na.rm = TRUE)
  n_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < -0.5, na.rm = TRUE)
  
  message(sprintf("  %s: %d sig (↑%d, ↓%d)", coef_name, n_sig, n_up, n_down))
  
  return(tt)
}

# Extract all contrasts
contrasts_list <- colnames(contrast_matrix)
results_list <- list()

for (contrast in contrasts_list) {
  results_list[[contrast]] <- save_limma_results(fit2, contrast)
}

# --------------------------------------------------------
# STEP 12: Volcano Plots
# --------------------------------------------------------
message("Creating volcano plots...")

save_volcano <- function(tt, contrast_name, top_n = 15) {
  tt <- tt %>%
    mutate(
      sig = case_when(
        adj.P.Val < 0.05 & logFC > 0.5 ~ "Up",
        adj.P.Val < 0.05 & logFC < -0.5 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Top genes to label
  top_up <- tt %>% filter(sig == "Up") %>% arrange(adj.P.Val) %>% head(top_n/2)
  top_down <- tt %>% filter(sig == "Down") %>% arrange(adj.P.Val) %>% head(top_n/2)
  top_labels <- bind_rows(top_up, top_down)
  
  p <- ggplot(tt, aes(x = logFC, y = -log10(P.Value), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_labels, 
                    aes(label = protein_id), 
                    size = 2.5,
                    max.overlaps = 20) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey70")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    theme_minimal(base_size = 12) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      title = paste("Volcano:", gsub("_", " ", contrast_name))
    )
  
  ggsave(file.path(OUT_DIR, "figs", paste0("volcano_", contrast_name, ".png")),
         p, width = 10, height = 8, dpi = 200)
}

# Create volcano for key comparisons
save_volcano(results_list[["TauWT_vs_Control_Day01"]], "TauWT_vs_Control_Day01")
save_volcano(results_list[["TauWT_vs_Control_Day10"]], "TauWT_vs_Control_Day10")
save_volcano(results_list[["TauWT_vs_Control_Day20"]], "TauWT_vs_Control_Day20")
save_volcano(results_list[["TauR406W_vs_Control_Day20"]], "TauR406W_vs_Control_Day20")

message("✓ Volcano plots saved")

# --------------------------------------------------------
# STEP 13: Heatmaps for All Key Contrasts
# --------------------------------------------------------
message("Creating heatmaps for all contrasts...")

create_heatmap <- function(contrast_name, results_df, expr_mat, prot_meta, 
                           top_n = 50, padj_cut = 0.1) {
  
  # Get top DE proteins
  top_proteins <- results_df %>%
    dplyr::filter(adj.P.Val < padj_cut) %>%
    arrange(adj.P.Val) %>%
    head(top_n) %>%
    pull(protein_id)
  
  if (length(top_proteins) < 5) {
    # If not enough significant, take top by p-value
    top_proteins <- results_df %>%
      arrange(P.Value) %>%
      head(top_n) %>%
      pull(protein_id)
    message("  ", contrast_name, ": Using top ", length(top_proteins), " by p-value (few significant)")
  } else {
    message("  ", contrast_name, ": ", length(top_proteins), " significant proteins")
  }
  
  # Subset expression matrix
  heat_mat <- expr_mat[top_proteins, ]
  
  # Z-score normalize
  heat_mat_z <- t(scale(t(heat_mat)))
  heat_mat_z[heat_mat_z > 2.5] <- 2.5
  heat_mat_z[heat_mat_z < -2.5] <- -2.5
  
  # Create simplified column labels (Genotype_Day_Rep)
  # Create replicate numbers within each group
  prot_meta_ordered <- prot_meta %>%
    group_by(genotype_simple, day) %>%
    mutate(rep = row_number()) %>%
    ungroup()
  
  col_labels <- paste0(
    substr(prot_meta_ordered$genotype_simple, 1, 4),  # Ctrl, TauW, TauR
    "_D", prot_meta_ordered$ageDays,
    "_", prot_meta_ordered$rep
  )
  
  # Annotation
  ann_col <- data.frame(
    Genotype = prot_meta$genotype_simple,
    Day = prot_meta$day,
    row.names = prot_meta$sample
  )
  
  ann_colors <- list(
    Genotype = c(Control = "steelblue", TauWT = "salmon", TauR406W = "darkred"),
    Day = c(Day01 = "#fee8c8", Day10 = "#fdbb84", Day20 = "#e34a33")
  )
  
  # Order columns by genotype then day
  col_order <- order(prot_meta$genotype_simple, prot_meta$day)
  heat_mat_z <- heat_mat_z[, col_order]
  ann_col <- ann_col[col_order, ]
  col_labels <- col_labels[col_order]
  
  # Simplify row labels (protein names)
  row_labels <- sapply(rownames(heat_mat_z), function(x) {
    # Extract gene name from protein ID (e.g., sp|P02515|HSP22_DROME -> HSP22)
    parts <- strsplit(x, "\\|")[[1]]
    if (length(parts) >= 3) {
      gene <- gsub("_DROME.*", "", parts[3])
      return(gene)
    }
    # If format different, truncate
    if (nchar(x) > 20) return(paste0(substr(x, 1, 17), "..."))
    return(x)
  })
  
  # Save heatmap
  png(file.path(OUT_DIR, "figs", paste0("heatmap_", contrast_name, ".png")),
      width = 1400, height = 1200, res = 150)
  
  pheatmap(heat_mat_z,
           annotation_col = ann_col,
           annotation_colors = ann_colors,
           labels_col = col_labels,
           labels_row = row_labels,
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 7,
           fontsize_col = 8,
           border_color = NA,
           main = paste0("Top DE Proteins: ", gsub("_", " ", contrast_name)))
  
  dev.off()
  message("  ✓ Saved: heatmap_", contrast_name, ".png")
}

# Create heatmaps for all key contrasts
key_contrasts <- c(
  "TauWT_vs_Control_Day01",
  "TauWT_vs_Control_Day10", 
  "TauWT_vs_Control_Day20",
  "TauR406W_vs_Control_Day01",
  "TauR406W_vs_Control_Day10",
  "TauR406W_vs_Control_Day20",
  "TauR406W_vs_TauWT_Day20"
)

for (contrast in key_contrasts) {
  if (contrast %in% names(results_list)) {
    create_heatmap(
      contrast_name = contrast,
      results_df = results_list[[contrast]],
      expr_mat = expr_imputed,
      prot_meta = prot_meta
    )
  }
}

message("✓ All heatmaps created!")

# --------------------------------------------------------
# STEP 14: Summary
# --------------------------------------------------------
message("\n========================================")
message("PROTEOMICS ANALYSIS COMPLETE!")
message("========================================")
message("\nOutput in: ", OUT_DIR)
message("\nDE Summary:")
for (contrast in contrasts_list) {
  n_sig <- sum(results_list[[contrast]]$adj.P.Val < 0.05, na.rm = TRUE)
  message(sprintf("  %s: %d significant proteins", contrast, n_sig))
}
message("\nFiles generated:")
message("  - tables/proteomics_metadata.csv")
message("  - tables/DE_*.csv (all contrasts)")
message("  - figs/01_PCA_proteomics.png")
message("  - figs/volcano_*.png")
message("  - figs/heatmap_top50_DE_proteins.png")

# ============================================================
# STEP 15: GSEA Analysis for Proteomics
# ============================================================

message("\n========================================")
message("Starting GSEA Analysis...")
message("========================================")

# --- Load additional packages ---
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")
if (!requireNamespace("org.Dm.eg.db", quietly = TRUE)) BiocManager::install("org.Dm.eg.db")
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db")

library(fgsea)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(GO.db)
library(data.table)

# --------------------------------------------------------
# --------------------------------------------------------
# STEP 16: Map Protein IDs to Gene Symbols (FIXED)
# --------------------------------------------------------
message("Mapping protein IDs to FlyBase gene symbols...")

# Extract UniProt accession from protein ID
# Format: sp|P02515|HSP22_DROME -> P02515
# Format: tr|Q9VZR2|Q9VZR2_DROME -> Q9VZR2

extract_uniprot_id <- function(protein_id) {
  parts <- strsplit(protein_id, "\\|")[[1]]
  if (length(parts) >= 2) {
    return(parts[2])
  }
  return(NA)
}

# Create mapping with UniProt IDs
protein_ids <- rownames(expr_imputed)

uniprot_map <- data.frame(
  protein_id = protein_ids,
  uniprot_id = sapply(protein_ids, extract_uniprot_id),
  stringsAsFactors = FALSE
)

message("  Extracted ", sum(!is.na(uniprot_map$uniprot_id)), " UniProt accessions")

# Get unique UniProt IDs
uniprot_ids <- unique(na.omit(uniprot_map$uniprot_id))

# Map to SYMBOL using org.Dm.eg.db
uniprot_to_symbol <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = uniprot_ids,
  keytype = "UNIPROT",
  columns = c("UNIPROT", "SYMBOL")
)

message("  Successfully mapped: ", sum(!is.na(uniprot_to_symbol$SYMBOL)), " proteins to FlyBase symbols")

# Merge back with protein IDs
gene_map <- uniprot_map %>%
  left_join(uniprot_to_symbol, by = c("uniprot_id" = "UNIPROT")) %>%
  dplyr::rename(gene_symbol = SYMBOL) %>%
  dplyr::filter(!is.na(gene_symbol) & gene_symbol != "") %>%
  distinct(protein_id, .keep_all = TRUE)

message("  Final gene map: ", nrow(gene_map), " proteins with FlyBase symbols")

# Check overlap with GO pathways
genes_in_pathways <- unique(unlist(pathways))
overlap <- sum(gene_map$gene_symbol %in% genes_in_pathways)
message("  Genes overlapping with GO pathways: ", overlap)

# --------------------------------------------------------
# STEP 17: Build Ranked Gene Lists for Each Contrast
# --------------------------------------------------------
message("Building ranked gene lists...")

build_ranks <- function(de_results, gene_map) {
  # Check if protein_id is already a column
  if ("protein_id" %in% colnames(de_results)) {
    de_with_genes <- de_results %>%
      left_join(gene_map, by = "protein_id")
  } else {
    de_with_genes <- de_results %>%
      rownames_to_column("protein_id") %>%
      left_join(gene_map, by = "protein_id")
  }
  
  # Filter valid entries
  de_with_genes <- de_with_genes %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::filter(!is.na(logFC), !is.na(P.Value))
  
  # For duplicate genes, keep the one with lowest p-value
  de_unique <- de_with_genes %>%
    group_by(gene_symbol) %>%
    slice_min(P.Value, n = 1) %>%
    ungroup()
  
  # Create rank score: sign(logFC) * -log10(pvalue)
  de_unique <- de_unique %>%
    mutate(rank_score = sign(logFC) * -log10(P.Value))
  
  # Named vector for fgsea
  ranks <- setNames(de_unique$rank_score, de_unique$gene_symbol)
  ranks <- sort(ranks, decreasing = TRUE)
  ranks <- ranks[is.finite(ranks)]
  
  return(ranks)
}

# --------------------------------------------------------
# STEP 18: Build GO Gene Sets
# --------------------------------------------------------
message("Building GO Biological Process gene sets...")

# Get all GO BP terms for Drosophila
go_bp <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = keys(org.Dm.eg.db, keytype = "SYMBOL"),
  columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"),
  keytype = "SYMBOL"
) %>%
  filter(ONTOLOGYALL == "BP") %>%
  dplyr::select(SYMBOL, GOALL) %>%
  distinct()

# Convert to list format for fgsea
pathways <- split(go_bp$SYMBOL, go_bp$GOALL)

# Filter by size (10-500 genes per pathway)
pathways <- pathways[sapply(pathways, length) >= 10 & sapply(pathways, length) <= 500]
message("GO BP pathways loaded: ", length(pathways))

# Get GO term names
go_terms <- AnnotationDbi::select(GO.db, keys = names(pathways), 
                                  columns = "TERM", keytype = "GOID")
term_map <- setNames(go_terms$TERM, go_terms$GOID)

# --------------------------------------------------------
# STEP 19: Run fGSEA for Each Contrast
# --------------------------------------------------------
run_gsea_for_contrast <- function(contrast_name, de_results, gene_map, pathways, term_map) {
  
  message("\n--- Running GSEA for: ", contrast_name, " ---")
  
  # Build ranks using the updated function
  ranks <- build_ranks(de_results, gene_map)
  message("  Genes in ranked list: ", length(ranks))
  
  if (length(ranks) < 50) {
    message("  WARNING: Too few genes for GSEA")
    return(NULL)
  }
  
  # Run fgsea
  set.seed(42)
  fgsea_res <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    nPermSimple = 10000
  )
  
  # Add GO term names
  fgsea_res <- fgsea_res %>%
    mutate(term = term_map[pathway]) %>%
    arrange(padj)
  
  n_sig <- sum(fgsea_res$padj < 0.05, na.rm = TRUE)
  message("  Significant pathways (padj < 0.05): ", n_sig)
  
  # Save results table
  fgsea_out <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
    dplyr::select(pathway, term, size, NES, pval, padj, leadingEdge)
  
  write.csv(fgsea_out, 
            file.path(GSEA_DIR, "tables", paste0("fgsea_GO_BP_", contrast_name, ".csv")),
            row.names = FALSE)
  
  # --- Barplot of top pathways ---
  top_paths <- fgsea_res %>%
    dplyr::filter(is.finite(padj)) %>%
    head(20) %>%
    mutate(label = ifelse(is.na(term) | term == "", pathway, term))
  
  if (nrow(top_paths) > 0) {
    p <- ggplot(top_paths, aes(x = reorder(label, NES), y = NES, fill = -log10(padj))) +
      geom_col() +
      coord_flip() +
      scale_fill_gradient(low = "steelblue", high = "darkred") +
      labs(
        title = paste0("GSEA GO Biological Process\n", gsub("_", " ", contrast_name)),
        x = NULL,
        y = "Normalized Enrichment Score (NES)",
        fill = "-log10(FDR)"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(file.path(GSEA_DIR, "figs", paste0("fgsea_GO_BP_", contrast_name, ".png")),
           p, width = 10, height = 8, dpi = 300)
  }
  
  # --- Mito-specific pathways ---
  mito_pat <- paste(
    "mitochond", "oxidative phosph", "electron transport",
    "respiratory chain", "ATP synthase", "tricarboxylic",
    "TCA", "beta-oxid", "mitophagy", "fission", "fusion",
    sep = "|"
  )
  
  mito_paths <- fgsea_res %>%
    dplyr::filter(grepl(mito_pat, term, ignore.case = TRUE)) %>%
    arrange(padj)
  
  if (nrow(mito_paths) > 0) {
    write.csv(
      mito_paths %>% mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")),
      file.path(GSEA_DIR, "tables", paste0("fgsea_mito_", contrast_name, ".csv")),
      row.names = FALSE
    )
    
    message("  Mito-related pathways: ", nrow(mito_paths))
    
    # Mito barplot
    mito_plot <- mito_paths %>%
      head(15) %>%
      mutate(label = ifelse(is.na(term) | term == "", pathway, term))
    
    if (nrow(mito_plot) > 0) {
      p_mito <- ggplot(mito_plot, aes(x = reorder(label, NES), y = NES, fill = -log10(padj))) +
        geom_col() +
        coord_flip() +
        scale_fill_gradient(low = "steelblue", high = "darkred") +
        labs(
          title = paste0("Mitochondrial Pathways\n", gsub("_", " ", contrast_name)),
          x = NULL,
          y = "NES",
          fill = "-log10(FDR)"
        ) +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(file.path(GSEA_DIR, "figs", paste0("fgsea_mito_", contrast_name, ".png")),
             p_mito, width = 10, height = 6, dpi = 300)
    }
  } else {
    message("  No mito-related pathways found")
  }
  
  return(fgsea_res)
}

# --------------------------------------------------------
# STEP 20: Run GSEA for All Contrasts
# --------------------------------------------------------
message("\nRunning GSEA for all contrasts...")

gsea_results <- list()

for (contrast in names(results_list)) {
  gsea_results[[contrast]] <- run_gsea_for_contrast(
    contrast_name = contrast,
    de_results = results_list[[contrast]],
    gene_map = gene_map,
    pathways = pathways,
    term_map = term_map
  )
}

# --------------------------------------------------------
# STEP 21: Summary
# --------------------------------------------------------
message("\n========================================")
message("GSEA ANALYSIS COMPLETE!")
message("========================================")
message("\nOutput in: ", normalizePath(GSEA_DIR))
message("\nGSEA Summary:")
for (contrast in names(gsea_results)) {
  if (!is.null(gsea_results[[contrast]])) {
    n_sig <- sum(gsea_results[[contrast]]$padj < 0.05, na.rm = TRUE)
    n_mito <- sum(grepl("mitochond|oxidative|electron|respiratory|ATP|TCA", 
                        gsea_results[[contrast]]$term, ignore.case = TRUE) & 
                    gsea_results[[contrast]]$padj < 0.05, na.rm = TRUE)
    message(sprintf("  %s: %d significant pathways (%d mito-related)", 
                    contrast, n_sig, n_mito))
  }
}
message("\nFiles generated:")
message("  - tables/fgsea_GO_BP_*.csv (all pathways)")
message("  - tables/fgsea_mito_*.csv (mito pathways)")
message("  - figs/fgsea_GO_BP_*.png (top pathway barplots)")
message("  - figs/fgsea_mito_*.png (mito pathway barplots)")

# --- Create GSEA trend plots ---
message("\nCreating GSEA trend plots...")

for (contrast in names(gsea_results)) {
  if (!is.null(gsea_results[[contrast]]) && nrow(gsea_results[[contrast]]) > 0) {
    
    # Convert to data.frame explicitly
    fg <- as.data.frame(gsea_results[[contrast]])
    
    # Get top 10 positive and top 10 negative NES
    top_up <- fg[fg$NES > 0, ]
    top_up <- top_up[order(top_up$pval), ]
    top_up <- head(top_up, 10)
    
    top_down <- fg[fg$NES < 0, ]
    top_down <- top_down[order(top_down$pval), ]
    top_down <- head(top_down, 10)
    
    top_paths <- rbind(top_up, top_down)
    
    if (nrow(top_paths) > 0) {
      # Create short term labels
      top_paths$term_short <- ifelse(
        is.na(top_paths$term) | top_paths$term == "", 
        top_paths$pathway,
        ifelse(nchar(top_paths$term) > 45, 
               paste0(substr(top_paths$term, 1, 42), "..."), 
               top_paths$term)
      )
      
      p <- ggplot(top_paths, 
                  aes(x = reorder(term_short, NES), y = NES, fill = -log10(padj))) +
        geom_col() +
        coord_flip() +
        scale_fill_gradient(low = "steelblue", high = "darkred", 
                            name = "-log10(FDR)") +
        theme_minimal(base_size = 10) +
        labs(
          x = "",
          y = "Normalized Enrichment Score (NES)",
          title = paste0("GSEA: ", gsub("_", " ", contrast)),
          subtitle = "Top pathways | Positive NES = Upregulated"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        theme(axis.text.y = element_text(size = 8))
      
      out_file <- file.path(GSEA_DIR, "figs", paste0("fgsea_trends_", contrast, ".png"))
      ggsave(out_file, p, width = 11, height = 8, dpi = 200)
    }
  }
}

message("✓ GSEA trend plots created")