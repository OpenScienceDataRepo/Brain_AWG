# Procedure

## QC and DESeq2

The raw counts were not available. When I reviewed the original paper, however, I found the following supplementary file:

*Table S1. Differentially expressed genes identified from different comparisons under FDR ≤ 0.05 and absolute log2 fold change ≥ 1.2. (XLSX 1236 kb)*

This meant that the QC and DESeq2 was already completed by the original researchers, making the data ready for the rest of the workflow.

*Sekiya M, Wang M, Fujisaki N, Sakakibara Y, Quan X, Ehrlich ME, De Jager PL, Bennett DA, Schadt EE, Gandy S, Ando K, Zhang B, Iijima KM. Integrated biology approach reveals molecular and pathological interactions among Alzheimer's Aβ42, Tau, TREM2, and TYROBP in Drosophila models. Genome Med. 2018 Mar 29;10(1):26. doi: 10.1186/s13073-018-0530-9. PMID: 29598827; PMCID: PMC5875009.*

## Visualization

To focus specifically on gene expression changes relative to baseline conditions, only contrasts that directly compared experimental groups to the control were selected. These were:
1. eGRL_TREM2-WT/TYROBP vs eGRL_Control
2. eGRL_TREM2-R47H/TYROBP vs eGRL_Control
3. eGRL_Aß42 vs eGRL_Control
4. eGRL_Aß42/TREM2-WT/TYROBP vs eGRL_Control
5.eGRL_Aß42/TREM2-R47H/TYROBP vs eGRL_Control

These contrasts were prioritized because they isolate the transcriptional effects of each genetic modification or pathological stimulus relative to the unaltered baseline, allowing clearer interpretation of up- or downregulated pathways under each condition.

All DEGs were plotted for those contrasts, and only mitochondrially relecant genes were labelled for interpretation.

Note: We consolidated all five contrasts into a single heatmap of mitochondrial genes. Attempting separate heatmaps for each contrast often produced too few genes for meaningful clustering. The combined heatmap provides a clear overview of expression patterns across all conditions, highlighting genes consistently up- or down-regulated while ensuring robust and interpretable clustering.

```
# Setup Working Directory
setwd("~/Desktop/ADBR Mito/GSE99012")  # change if needed

# Load Packages
pkgs_cran <- c("ggplot2","ggrepel","pheatmap","dplyr","readxl","tidyr")
for (p in pkgs_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("AnnotationDbi","org.Dm.eg.db","GO.db")  # fly annotation
for (p in pkgs_bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(readxl)
  library(tidyr)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  library(GO.db)
})

# Output Directory
OUT_DIR <- "RESULTS_GSE99012"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Load DEG Table (Table S1)
deg_table <- read_excel("DEG Sekiya et al. (2018).xlsx", sheet = "Table S1")
deg_table <- deg_table[-c(1,2), ]  # Remove metadata/header rows
deg_table <- as.data.frame(deg_table)  # ensure data frame

# Rename columns
colnames(deg_table)[1:8] <- c("GeneSymbol","Symbol","human_ortholog","log2FC","t","PValue","FDR","Comparison")

# Convert numeric columns
deg_table$log2FC <- as.numeric(deg_table$log2FC)
deg_table$FDR    <- as.numeric(deg_table$FDR)
deg_table$GeneSymbol <- as.character(deg_table$GeneSymbol)
deg_table$Symbol <- as.character(deg_table$Symbol)
deg_table$Comparison <- as.character(deg_table$Comparison)

# Curated Drosophila Mitochondrial Genes (FBgn IDs via GO)
mito_go_bp <- c("GO:0006119","GO:0022900","GO:0006099","GO:0006635","GO:0000422")
mito_go_cc <- c("GO:0005739","GO:0005743","GO:0005747","GO:0005753","GO:0005759")

.expand_offspring <- function(go_ids, offspr_env) {
  kids <- unique(unlist(mget(go_ids, envir=offspr_env, ifnotfound=NA)))
  unique(na.omit(c(go_ids, kids)))
}

all_go <- unique(c(.expand_offspring(mito_go_bp, GOBPOFFSPRING),
                   .expand_offspring(mito_go_cc, GOCCOFFSPRING)))

# Map GO → FBgn IDs
go2genes <- AnnotationDbi::select(org.Dm.eg.db, keys=all_go, keytype="GO", columns="FLYBASE")
mito_genes <- unique(na.omit(go2genes$FLYBASE))

# Volcano Plot Function
save_volcano_mito <- function(df, tag, mito_genes, thr_FDR=0.05, thr_LFC=1) {
  
  if(all(is.na(df$log2FC)) | all(is.na(df$FDR))) {
    message("Skipping ", tag, ": no numeric log2FC/FDR")
    return(NULL)
  }
  
  df$group <- "NS"
  df$group[df$FDR < thr_FDR & df$log2FC >=  thr_LFC] <- "Up (sig)"
  df$group[df$FDR < thr_FDR & df$log2FC <= -thr_LFC] <- "Down (sig)"
  df$group <- factor(df$group, levels=c("Down (sig)","NS","Up (sig)"))
  
  lab_set <- df$GeneSymbol[df$GeneSymbol %in% mito_genes & df$group != "NS"]
  
  out <- file.path(OUT_DIR, paste0("volcano_", gsub("[ /]","_",tag), ".png"))
  png(out, width=1800, height=1500, res=200, bg="white")
  
  print(
    ggplot(df, aes(x=log2FC, y=-log10(FDR), color=group)) +
      geom_point(alpha=0.7, size=1.6, na.rm=TRUE) +
      geom_point(data=subset(df, GeneSymbol %in% lab_set), shape=21, stroke=0.7, size=2.8, color="black") +
      geom_vline(xintercept=c(-thr_LFC, thr_LFC), linetype="dashed") +
      geom_hline(yintercept=-log10(thr_FDR), linetype="dashed") +
      ggrepel::geom_text_repel(
        data=subset(df, GeneSymbol %in% lab_set),
        aes(label=Symbol),
        size=3, max.overlaps=500, box.padding=0.4, min.segment.length=0
      ) +
      scale_color_manual(values=c("Down (sig)"="#2C7BB6","NS"="grey60","Up (sig)"="#D7191C")) +
      labs(x="log2 Fold Change", y="-log10(FDR)",
           title=paste0("Volcano (mitochondrial focus): ", tag),
           subtitle="Dashed: |log2FC|=1 & FDR=0.05; Labels = significant mito genes only.") +
      theme_classic() + theme(legend.position="top")
  )
  
  dev.off()
}

# Key Comparisons: Only conditions vs Control
comparisons_to_plot <- c(
  "eGRL_TREM2-WT/TYROBP vs eGRL_Control",
  "eGRL_TREM2-R47H/TYROBP vs eGRL_Control",
  "eGRL_Aß42 vs eGRL_Control",
  "eGRL_Aß42/TREM2-WT/TYROBP vs eGRL_Control",
  "eGRL_Aß42/TREM2-R47H/TYROBP vs eGRL_Control"
)

# Generate Volcano Plots for Each Contrast
for (comp in comparisons_to_plot) {
  df_sub <- deg_table[deg_table$Comparison == comp, ]
  save_volcano_mito(df_sub, tag=comp, mito_genes=mito_genes)
}

# Unified Heatmap Function (mitochondrial only, with row-check)
save_heatmap_unified <- function(tag, 
                                 data_source,     
                                 comparisons=NULL, 
                                 res_unshr=NULL,   
                                 meta=NULL,        
                                 topN=40, cap_z=2.5) {
  
  fbgn2symbol <- setNames(deg_table$Symbol, deg_table$GeneSymbol)
  
  if (data_source == "deg_table") {
    if (is.null(comparisons)) stop("comparisons must be provided for deg_table heatmap")
    
    # Filter to mitochondrial genes only
    mat_fc <- deg_table %>%
      dplyr::filter(
    Comparison %in% comparisons,
    GeneSymbol %in% mito_genes,
    FDR < 0.05,
    abs(log2FC) >= 1
  ) %>%
  dplyr::select(GeneSymbol, log2FC, Comparison) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = log2FC)
    
    # Fill NAs with 0 for plotting
mat_fc[is.na(mat_fc)] <- 0

# Rank by max |log2FC|
mat_fc$max_abs_LFC <- apply(abs(mat_fc[,-1]), 1, max, na.rm=TRUE)

mat_fc_top <- mat_fc %>%
  arrange(desc(max_abs_LFC)) %>%
  slice_head(n = min(topN, nrow(mat_fc)))

# Row labels
row_labels <- fbgn2symbol[mat_fc_top$GeneSymbol]
row_labels[is.na(row_labels)] <- mat_fc_top$GeneSymbol

# Plotting matrix (exclude GeneSymbol and max_abs_LFC)
plot_mat <- as.matrix(mat_fc_top[, !(colnames(mat_fc_top) %in% c("GeneSymbol","max_abs_LFC"))])

    
  } else if (data_source == "DESeq2") {
    if (is.null(res_unshr) | is.null(meta)) stop("res_unshr and meta must be provided for DESeq2 heatmap")
    
    mat <- assay(vsd)
    rdf <- as.data.frame(res_unshr)
    rdf$gene_id <- rownames(rdf)
    rdf <- rdf[order(rdf$padj, -abs(rdf$log2FoldChange)), ]
    genes_use <- head(stats::na.omit(rdf$gene_id), topN)
    genes_use <- intersect(genes_use, rownames(mat))
    if (!length(genes_use)) stop("No genes available for heatmap: ", tag)
    
    cols_use <- rownames(meta)
    plot_mat <- mat[genes_use, cols_use, drop=FALSE]
    
    vars <- apply(plot_mat, 1, var)
    plot_mat <- plot_mat[vars > 0, , drop=FALSE]
    
    if(nrow(plot_mat) < 2){
      message("Skipping heatmap ", tag, ": less than 2 genes with variable expression")
      return(NULL)
    }
    
    row_labels <- fbgn2symbol[rownames(plot_mat)]
    row_labels[is.na(row_labels)] <- rownames(plot_mat)
    
  } else {
    stop("Invalid data_source: choose 'deg_table' or 'DESeq2'")
  }
  
  # Row-Z scale
  z <- t(scale(t(plot_mat)))
  z[is.na(z)] <- 0
  z[z > cap_z] <- cap_z
  z[z < -cap_z] <- -cap_z
  
  ann <- NULL
  if (data_source == "DESeq2") {
    ann <- meta[, c("condition_group","sex"), drop=FALSE]
    ann$replicate <- if("replicate" %in% colnames(meta)) meta$replicate else seq_len(nrow(ann))
  }
  
  png(file.path(OUT_DIR, paste0("heatmap_", gsub("[ /]","_",tag), ".png")), 
      width=1800, height=1500, res=200, bg="white")
  pheatmap(z,
           labels_row = row_labels,
           show_rownames = TRUE,
           cluster_rows=TRUE,
           cluster_cols=TRUE,
           annotation_col=ann,
           main=paste0("Top ", nrow(z), " mito genes (row-Z cap ±", cap_z, ")"),
           fontsize_row=7, fontsize_col=9,
           border_color=NA)
  dev.off()
}

# Generate Unified Heatmap Across All Contrasts
save_heatmap_unified("log2FC_vsControl_mito", data_source="deg_table", comparisons=comparisons_to_plot)

message("✅ Saved ", OUT_DIR)

```
# GSEA

```
# Set working directory (change if needed)
setwd("~/Desktop/ADBR Mito/GSE99012")

# Packages
pkgs_cran <- c("data.table","ggplot2")
for (p in pkgs_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("fgsea","AnnotationDbi","org.Dm.eg.db","GO.db")
for (p in pkgs_bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(fgsea);      library(AnnotationDbi)
  library(org.Dm.eg.db); library(GO.db)
})

# Output dirs
GSEA_DIR <- "GSEA"
dir.create(GSEA_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(GSEA_DIR,"tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(GSEA_DIR,"figs"),   recursive=TRUE, showWarnings=FALSE)

# Load DEG table (Sekiya et al. Table S1)
deg_table <- readxl::read_excel("DEG Sekiya et al. (2018).xlsx", sheet = "Table S1")
deg_table <- deg_table[-c(1,2), ]
deg_table <- as.data.frame(deg_table)
colnames(deg_table)[1:8] <- c("GeneSymbol","Symbol","human_ortholog","log2FC","t","PValue","FDR","Comparison")

# === FIX: force numeric conversion ===
deg_table$log2FC <- as.numeric(gsub(",", ".", deg_table$log2FC))
deg_table$PValue <- as.numeric(gsub(",", ".", deg_table$PValue))
deg_table$FDR    <- as.numeric(gsub(",", ".", deg_table$FDR))

deg_table$GeneSymbol <- as.character(deg_table$GeneSymbol)
deg_table$Symbol     <- as.character(deg_table$Symbol)
deg_table$Comparison <- as.character(deg_table$Comparison)

# List of contrasts
comparisons_to_plot <- c(
  "eGRL_TREM2-WT/TYROBP vs eGRL_Control",
  "eGRL_TREM2-R47H/TYROBP vs eGRL_Control",
  "eGRL_Aß42 vs eGRL_Control",
  "eGRL_Aß42/TREM2-WT/TYROBP vs eGRL_Control",
  "eGRL_Aß42/TREM2-R47H/TYROBP vs eGRL_Control"
)

# --- Helpers (unchanged) ---
build_ranks <- function(res_df) {
  if (!"GeneSymbol" %in% colnames(res_df)) stop("res_df must contain GeneSymbol column")
  ids <- res_df$GeneSymbol
  if ("stat" %in% names(res_df) && any(is.finite(res_df$stat))) {
    score <- res_df$stat
  } else {
    pv <- if ("PValue" %in% names(res_df)) res_df$PValue else res_df$FDR
    pv[!is.finite(pv) | is.na(pv)] <- 1
    lfc <- if ("log2FC" %in% names(res_df)) res_df$log2FC else 0
    score <- sign(lfc) * (-log10(pmax(pv, 1e-300)))
  }
  keep <- is.finite(score) & !is.na(ids) & nzchar(ids)
  ranks <- tapply(score[keep], ids[keep], mean)
  sort(ranks, decreasing=TRUE)
}

build_go_bp_sets <- function(ranks, min_size=10, max_size=500) {
  genes <- names(ranks)
  m <- AnnotationDbi::select(org.Dm.eg.db,
                             keys=genes, keytype="FLYBASE",
                             columns=c("GOALL","ONTOLOGYALL"))
  m <- unique(as.data.frame(m, stringsAsFactors=FALSE))
  m <- m[!is.na(m$GOALL) & !is.na(m$ONTOLOGYALL) & m$ONTOLOGYALL=="BP", , drop=FALSE]
  pathways <- split(m$FLYBASE, m$GOALL, drop=TRUE)
  pathways <- lapply(pathways, function(x) unique(x[!is.na(x) & nzchar(x)]))
  lens <- vapply(pathways, length, integer(1))
  pathways[lens >= min_size & lens <= max_size]
}

term_map_for <- function(go_ids) {
  if (!length(go_ids)) return(setNames(character(0), character(0)))
  tbl <- AnnotationDbi::select(GO.db, keys=go_ids, keytype="GOID", columns="TERM")
  tbl <- unique(as.data.frame(tbl, stringsAsFactors=FALSE))
  setNames(tbl$TERM, tbl$GOID)
}

`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

pick_mito_paths <- function(fg_dt, padj_cut=0.25, top_n=8) {
  if (!nrow(fg_dt)) return(character(0))
  lab <- if ("term" %in% names(fg_dt)) fg_dt$term else fg_dt$pathway
  mito_pat <- paste(
    "mitochond", "oxidative phosph", "electron transport", "respiratory chain",
    "ATP synthase", "tricarboxylic", "TCA", "beta-oxid", "mitophagy", "fission", "fusion",
    sep="|"
  )
  mito <- fg_dt[grepl(mito_pat, lab, ignore.case=TRUE) & is.finite(fg_dt$padj), , drop=FALSE]
  mito <- mito[order(mito$padj, -abs(mito$NES)), , drop=FALSE]
  sel <- mito$pathway[mito$padj <= padj_cut]
  if (!length(sel)) sel <- head(mito$pathway, min(top_n, nrow(mito)))
  unique(na.omit(sel))
}

save_enrichment_curve <- function(tag, pathway_id, ranks, pathways, TERM_MAP) {
  if (!pathway_id %in% names(pathways)) return(invisible(NULL))
  term_label <- TERM_MAP[[pathway_id]] %||% pathway_id
  gp <- fgsea::plotEnrichment(pathways[[pathway_id]], ranks) +
        ggplot2::ggtitle(paste0(term_label, " (", tag, ")"))
  fn <- file.path(GSEA_DIR,"figs",
                  paste0("fgsea_curve_",
                         gsub("[^A-Za-z0-9]+","_", term_label), "_", gsub("[^A-Za-z0-9]+","_", tag), ".png"))
  ggsave(fn, gp, width=7, height=5, dpi=300, bg="white")
  message("Saved: ", fn)
}

save_gsea_table_plot <- function(tag, sel_paths, ranks, pathways, fg_dt, TERM_MAP) {
  sel_paths <- intersect(sel_paths, names(pathways))
  if (!length(sel_paths)) { message("No selected pathways to plot for ", tag); return(invisible(NULL)) }
  fg_sel <- fg_dt[match(sel_paths, fg_dt$pathway), , drop=FALSE]
  fn <- file.path(GSEA_DIR, "figs", paste0("fgsea_table_", gsub("[^A-Za-z0-9]+","_", tag), ".png"))
  png(fn, width=1600, height=900, res=150, bg="white")
  fgsea::plotGseaTable(pathways[sel_paths], ranks, fg_sel)
  dev.off()
  message("Saved: ", fn)
}

run_fgsea_and_plots <- function(tag, nperm=10000, padj_cut_table=0.25, topN_bar=15) {
  message("\n=== ", tag, " ===")
  res_df <- deg_table[deg_table$Comparison == tag, , drop=FALSE]
  if (!nrow(res_df)) { message("No DE rows for ", tag); return(invisible(NULL)) }
  ranks <- build_ranks(res_df)
  if (length(ranks) < 50) { message("Too few ranked genes for ", tag, "; skipping"); return(invisible(NULL)) }

  pathways <- build_go_bp_sets(ranks)
  if (!length(pathways)) { message("No GO BP pathways built for ", tag); return(invisible(NULL)) }
  TERM_MAP <- term_map_for(names(pathways))

  set.seed(42)
  fg <- suppressWarnings(fgsea(pathways=pathways, stats=ranks, nperm=nperm))
  fg_dt <- as.data.table(fg)
  if (!nrow(fg_dt)) { message("fgsea returned no results for ", tag); return(invisible(NULL)) }
  setorder(fg_dt, padj, -NES)
  fg_dt[, term := TERM_MAP[pathway]]

  out_csv <- file.path(GSEA_DIR, "tables", paste0("fgsea_GO_BP_", gsub("[^A-Za-z0-9]+","_", tag), ".csv"))
  fwrite(fg_dt[, .(pathway, term, size, NES, pval, padj,
                   leadingEdge = vapply(leadingEdge, \(v) paste(v, collapse=";"), character(1)))],
         out_csv)
  message("Saved: ", out_csv)

  # Top GO BP barplot
  topN <- head(fg_dt[is.finite(padj)], topN_bar)
  if (nrow(topN)) {
    topN[, label := ifelse(is.na(term) | term=="", pathway, term)]
    p <- ggplot(topN, aes(x=reorder(label, NES), y=NES, fill=-log10(padj))) +
      geom_col() + coord_flip() +
      labs(title=paste0("fgsea GO BP: ", gsub("_"," ", tag)),
           x=NULL, y="NES", fill="-log10(FDR)") +
      theme_minimal(base_size=12)
    out_png <- file.path(GSEA_DIR, "figs", paste0("fgsea_GO_BP_", gsub("[^A-Za-z0-9]+","_", tag), ".png"))
    ggsave(out_png, p, width=8, height=6, dpi=300, bg="white")
    message("Saved: ", out_png)
  }

  sel <- pick_mito_paths(fg_dt, padj_cut=padj_cut_table, top_n=8)
  if (!length(sel)) {
    message("No mito keywords found at the chosen FDR; falling back to top 6 terms for curves.")
    sel <- head(fg_dt$pathway, 6)
  }

  invisible(lapply(sel, save_enrichment_curve,
                   tag=tag, ranks=ranks, pathways=pathways, TERM_MAP=TERM_MAP))

  save_gsea_table_plot(tag, sel, ranks, pathways, fg_dt, TERM_MAP)
}

# Run for all contrasts
for (tag in comparisons_to_plot) {
  run_fgsea_and_plots(tag, nperm=10000)
}

message("\nAll done. See: ", normalizePath(GSEA_DIR))
```
