# 0. Setup
setwd("/Volumes/Marians_SSD/ADBR_Mito/PRJNA747152")

# Load packages
pkgs_cran <- c("ggplot2","ggrepel","pheatmap","dplyr")
for (p in pkgs_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("DESeq2","AnnotationDbi","org.Dm.eg.db","GO.db","apeglm")
for (p in pkgs_bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(pheatmap); library(dplyr)
  library(DESeq2); library(AnnotationDbi); library(org.Dm.eg.db); library(GO.db); library(apeglm)
})

# Output directories
OUT_DIR <- "DE_output"
dir.create(file.path(OUT_DIR,"figs"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR,"tables"), recursive=TRUE, showWarnings=FALSE)

# 1. Load data
counts <- read.delim("counts_gene_length_clean.txt", row.names=1, check.names=FALSE)
counts_numeric <- counts[, 7:ncol(counts)]

meta <- read.csv("meta.csv", stringsAsFactors=FALSE)
meta <- meta[meta$Run %in% colnames(counts_numeric), ]
rownames(meta) <- meta$Run

# 2. Define condition groups
meta$condition_group <- sapply(meta$treatment, function(x) {
  if (x == "0 Gy gamma rays") return("Control")
  if (x == "0.4 Gy gamma rays") return("LowGamma")
  if (x == "10 Gy gamma rays") return("HighGamma")
  if (x == "0.4+10Gy gamma rays") return("ComboGamma")
  return(NA)
})
meta$condition_group <- factor(meta$condition_group,
                               levels=c("Control","LowGamma","HighGamma","ComboGamma"))

# 3. DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_numeric,
  colData   = meta,
  design    = ~ condition_group
)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 0, ]

# Run DESeq
dds <- DESeq(dds)

# 4. Contrast function
run_contrast <- function(var, levelA, levelB, tag) {
  res <- results(dds, contrast=c(var, levelA, levelB), alpha=0.05)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  write.csv(res_df, file.path(OUT_DIR,"tables",paste0("DE_",tag,"_unshrunk.csv")), row.names=FALSE)
  
  # Shrink LFC
  coef_name <- paste0(var, "_", levelA, "_vs_", levelB)
  res_shr_df <- NULL
  if (coef_name %in% resultsNames(dds)) {
    res_shr <- lfcShrink(dds, coef=coef_name, type="apeglm")
    res_shr_df <- as.data.frame(res_shr)
    res_shr_df$gene <- rownames(res_shr)
    write.csv(res_shr_df,
              file.path(OUT_DIR,"tables",paste0("DE_",tag,"_shrunk.csv")),
              row.names=FALSE)
  }
  
  # Summary
  sig <- subset(res_df, !is.na(padj) & padj < 0.05)
  sig_lfc <- subset(sig, abs(log2FoldChange) >= 1)
  message(sprintf("[%s] padj<0.05: %d | padj<0.05 & |LFC|>=1: %d",
                  tag, nrow(sig), nrow(sig_lfc)))
  
  return(list(res=res_df, res_shr=res_shr_df))
}

# 5. Run contrasts
out_Low_vs_Control   <- run_contrast("condition_group","LowGamma","Control","Low_vs_Control")
out_High_vs_Control  <- run_contrast("condition_group","HighGamma","Control","High_vs_Control")
out_Combo_vs_Control <- run_contrast("condition_group","ComboGamma","Control","Combo_vs_Control")
out_Combo_vs_High    <- run_contrast("condition_group","ComboGamma","HighGamma","Combo_vs_High")
out_Combo_vs_Low     <- run_contrast("condition_group","ComboGamma","LowGamma","Combo_vs_Low")

# 6. Mitochondrial gene set
mito_go_bp <- c("GO:0006119","GO:0022900","GO:0006099","GO:0006635","GO:0000422")
mito_go_cc <- c("GO:0005739","GO:0005743","GO:0005747","GO:0005753","GO:0005759")

.expand_offspring <- function(go_ids, offspr_env) {
  kids <- unique(unlist(mget(go_ids, envir=offspr_env, ifnotfound=NA)))
  unique(na.omit(c(go_ids, kids)))
}

all_go <- unique(c(
  .expand_offspring(mito_go_bp, GOBPOFFSPRING),
  .expand_offspring(mito_go_cc, GOCCOFFSPRING)
))

go2genes <- AnnotationDbi::select(org.Dm.eg.db,
                                  keys=all_go,
                                  keytype="GO",
                                  columns="FLYBASE")
mito_fbgn <- unique(na.omit(go2genes$FLYBASE))

# 7. ID → symbol
toSym <- function(x) {
  map <- AnnotationDbi::select(org.Dm.eg.db,
                               keys=x,
                               keytype="FLYBASE",
                               columns="SYMBOL")
  sym <- setNames(map$SYMBOL, map$FLYBASE)
  y <- sym[x]
  y[is.na(y)] <- x[is.na(y)]
  return(y)
}

# 8. Volcano function
save_volcano <- function(res_df, tag, mito_fbgn,
                         thr_p=0.05, thr_fc=1) {
  vdf <- res_df
  vdf <- vdf[!is.na(vdf$padj), ]
  vdf$padj[vdf$padj == 0] <- 1e-300
  vdf$gene_id <- vdf$gene
  vdf$gene_symbol <- toSym(vdf$gene_id)
  vdf$log2FC <- vdf$log2FoldChange
  vdf$group <- "NS"
  vdf$group[vdf$padj < thr_p & vdf$log2FC >= thr_fc] <- "Up"
  vdf$group[vdf$padj < thr_p & vdf$log2FC <= -thr_fc] <- "Down"
  
  lab_set <- intersect(mito_fbgn, vdf$gene_id)
  if(length(lab_set)==0) lab_set <- head(vdf[order(vdf$padj), "gene_id"], 12)
  
  out <- file.path(OUT_DIR,"figs",paste0("volcano_",tag,".png"))
  png(out, width=1800, height=1500, res=200)
  print(
    ggplot(vdf, aes(x=log2FC, y=-log10(padj), color=group)) +
      geom_point(alpha=0.7, size=1.5) +
      geom_vline(xintercept=c(-thr_fc,thr_fc), linetype="dashed") +
      geom_hline(yintercept=-log10(thr_p), linetype="dashed") +
      geom_point(data=subset(vdf, gene_id %in% lab_set),
                 shape=21, color="black", size=2.5, stroke=0.6) +
      ggrepel::geom_text_repel(data=subset(vdf, gene_id %in% lab_set),
                               aes(label=gene_symbol),
                               size=3, max.overlaps=40) +
      scale_color_manual(values=c("Down"="#2C7BB6","NS"="grey70","Up"="#D7191C")) +
      labs(title=paste("Volcano:", tag),
           x="Log2 Fold Change",
           y="-log10(FDR)") +
      theme_classic() +
      theme(legend.position="top")
  )
  dev.off()
}

# 9. Generate volcano plots
save_volcano(out_Low_vs_Control$res,   "Low_vs_Control", mito_fbgn)
save_volcano(out_High_vs_Control$res,  "High_vs_Control", mito_fbgn)
save_volcano(out_Combo_vs_Control$res, "Combo_vs_Control", mito_fbgn)
save_volcano(out_Combo_vs_High$res,    "Combo_vs_High", mito_fbgn)
save_volcano(out_Combo_vs_Low$res,     "Combo_vs_Low", mito_fbgn)

# 10. PCA and boxplots
vsd <- vst(dds, blind=TRUE)

# PCA
pca_plot <- plotPCA(vsd, intgroup="condition_group") + ggtitle("PCA: Condition Groups") + theme_classic()
png(file.path(OUT_DIR, "figs", "PCA_condition_group.png"), width=1600, height=1200, res=200)
print(pca_plot)
dev.off()

# Normalized counts boxplot
norm_counts <- counts(dds, normalized=TRUE)
png(file.path(OUT_DIR, "figs", "boxplot_normalized_counts.png"), width=1800, height=1200, res=200)
boxplot(log2(norm_counts + 1),
        las=2,
        col="lightblue",
        main="Normalized Counts (log2 scale)",
        ylab="log2(count + 1)")
dev.off()

# 11. Heatmaps
save_heatmap <- function(res_df, tag, top_n=50, mito_only=FALSE, mito_fbgn=NULL) {
  df <- res_df
  df <- df[!is.na(df$padj), ]
  
  if(mito_only && !is.null(mito_fbgn)) {
    df <- df[df$gene %in% mito_fbgn, ]
    if(nrow(df)==0) {
      warning(paste("No mitochondrial genes found for", tag))
      return(NULL)
    }
  }
  
  # Top N genes
  df_top <- df[order(df$padj), ][1:min(top_n, nrow(df)), ]
  
  # Normalized counts
  norm_mat <- counts(dds, normalized=TRUE)
  mat <- norm_mat[df_top$gene, , drop=FALSE]
  rownames(mat) <- toSym(rownames(mat))
  mat_scaled <- t(scale(t(mat)))
  
  # Order columns by condition_group
  col_order <- order(meta$condition_group)
  mat_scaled <- mat_scaled[, col_order]
  
  ann_col <- data.frame(
    Condition = meta$condition_group[col_order]
  )
  rownames(ann_col) <- rownames(meta)[col_order]
  
  # Condition colors
  ann_colors <- list(
    Condition = c(
      Control="#1f77b4", LowGamma="#ff7f0e",
      HighGamma="#2ca02c", ComboGamma="#d62728"
    )
  )
  
  out_file <- file.path(OUT_DIR,"figs",paste0("heatmap_", tag, ifelse(mito_only,"_mito",""), ".png"))
  
  png(out_file, width=2000, height=1600, res=200)
  pheatmap(mat_scaled,
           annotation_col = ann_col,
           annotation_colors = ann_colors,
           show_rownames = TRUE,
           fontsize_row = 6,
           cluster_rows = TRUE,        # cluster genes
           cluster_cols = FALSE,       # keep columns grouped by condition
           main = paste("Heatmap:", tag, ifelse(mito_only,"(Mitochondrial Genes)",""))
  )
  dev.off()
}

# 12. Generate heatmaps
contrast_list <- list(
  Low_vs_Control   = out_Low_vs_Control$res,
  High_vs_Control  = out_High_vs_Control$res,
  Combo_vs_Control = out_Combo_vs_Control$res,
  Combo_vs_High    = out_Combo_vs_High$res,
  Combo_vs_Low     = out_Combo_vs_Low$res
)

# Full heatmaps
for(tag in names(contrast_list)) {
  save_heatmap(contrast_list[[tag]], tag, top_n=50, mito_only=FALSE)
}

# Mitochondrial-specific heatmaps
for(tag in names(contrast_list)) {
  save_heatmap(contrast_list[[tag]], tag, top_n=50, mito_only=TRUE, mito_fbgn=mito_fbgn)
}

message("All heatmaps saved in ", file.path(OUT_DIR,"figs"))

## GSEA

# Packages

pkgs_cran <- c(
  "data.table","ggplot2","dplyr","ggrepel",
  "enrichR","stringr","tidyr","igraph","ggraph"
)

pkgs_bioc <- c(
  "fgsea","AnnotationDbi","org.Dm.eg.db","GO.db"
)

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly=TRUE)) {
    BiocManager::install(p, ask=FALSE, update=FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  library(GO.db)
  library(enrichR)
  library(stringr)
  library(tidyr)
  library(igraph)
  library(ggraph)
})

# PATHS

BASE_DIR <- "/Volumes/Marians_SSD/ADBR_Mito/PRJNA747152"

GSEA_DIR <- file.path(BASE_DIR, "GSEA")
TBL_DIR  <- file.path(BASE_DIR, "DE_output", "tables")

dir.create(file.path(GSEA_DIR, "tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(GSEA_DIR, "figs"), recursive=TRUE, showWarnings=FALSE)

# PARAMETERS

PADJ_CUT <- 0.05
LFC_CUT  <- 1

CONTRASTS <- c("Combo_vs_Control","Combo_vs_High","Combo_vs_Low", "High_vs_Control", "Low_vs_Control")

# HELPERS

`%||%` <- function(a,b) if (!is.null(a) && length(a)) a else b

wrap_terms <- function(x,w=40) stringr::str_wrap(x,w)

# LOAD DE

load_de <- function(tag) {
  f <- file.path(TBL_DIR, paste0("DE_", tag, "_unshrunk.csv"))
  if (!file.exists(f)) stop("Missing DE file: ", f)

  df <- read.csv(f, stringsAsFactors=FALSE, check.names=FALSE)

  if (!"gene" %in% colnames(df)) stop("No gene column")

  rownames(df) <- df$gene
  df
}

# RANKS

build_ranks <- function(res) {

  ids <- rownames(res)

  if ("stat" %in% names(res) && any(is.finite(res$stat))) {
    score <- res$stat
  } else {
    pv  <- res$pvalue %||% res$padj
    lfc <- res$log2FoldChange %||% 0

    pv[!is.finite(pv)]  <- 1
    lfc[!is.finite(lfc)] <- 0

    score <- sign(lfc) * (-log10(pmax(pv,1e-300)))
  }

  keep <- is.finite(score) & nzchar(ids)
  ranks <- tapply(score[keep], ids[keep], mean)

  sort(ranks, decreasing=TRUE)
}

# GO BP SETS (RESTORED)

build_go_bp_sets <- function(ranks) {

  genes <- names(ranks)

  m <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys=genes,
    keytype="FLYBASE",
    columns=c("GOALL","ONTOLOGYALL")
  )

  m <- m[m$ONTOLOGYALL=="BP", ]

  if (nrow(m)==0) return(NULL)

  paths <- split(m$FLYBASE, m$GOALL)
  paths <- lapply(paths, function(x) unique(na.omit(x)))

  lens <- lengths(paths)

  paths[lens >= 10 & lens <= 500]
}

# GO TERM MAP

get_go_terms <- function(ids) {

  if (!length(ids)) return(setNames(character(0),character(0)))

  tbl <- AnnotationDbi::select(
    GO.db,
    keys=ids,
    keytype="GOID",
    columns="TERM"
  )

  tbl <- unique(tbl)
  setNames(tbl$TERM, tbl$GOID)
}

# SIGNIFICANT GENES

get_sig <- function(res) {

  res <- res[
    !is.na(res$padj) &
    res$padj < PADJ_CUT &
    !is.na(res$log2FoldChange) &
    abs(res$log2FoldChange) >= LFC_CUT,
  ]

  unique(rownames(res))
}

# SYMBOL CONVERSION

fbgn_to_symbol <- function(ids) {

  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (!length(ids)) return(character(0))

  tbl <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys=ids,
    keytype="FLYBASE",
    columns="SYMBOL"
  )

  tbl <- tbl[!duplicated(tbl$FLYBASE), ]

  map <- setNames(tbl$SYMBOL, tbl$FLYBASE)

  out <- map[ids]
  out[is.na(out)] <- ids[is.na(out)]

  unique(out)
}

# FGSEA PLOT

save_fgsea <- function(dt, tag) {

  if (is.null(dt) || nrow(dt)==0) return()

  dt <- dt[order(padj)][1:min(15,nrow(dt))]

  term_map <- get_go_terms(dt$pathway)

  dt$term <- term_map[dt$pathway]
  dt$term[is.na(dt$term)] <- dt$pathway

  dt$label <- wrap_terms(dt$term)

  p <- ggplot(dt, aes(reorder(label, NES), NES, fill=-log10(padj))) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(title=paste("GSEA:",tag))

  ggsave(file.path(GSEA_DIR,"figs",paste0("fgsea_",tag,".png")),
         p,width=10,height=6)
}

save_enrichment_plots <- function(pathways, ranks, fg_dt, tag) {

  if (is.null(fg_dt) || nrow(fg_dt) == 0) return()

  # top pathways (by padj)
  top <- fg_dt[order(padj)][1:min(5, nrow(fg_dt))]

  term_map <- get_go_terms(top$pathway)

  for (i in seq_len(nrow(top))) {

    pw <- top$pathway[i]
    term_name <- term_map[pw]
    if (is.na(term_name)) term_name <- pw

    p <- plotEnrichment(pathways[[pw]], ranks) +
      labs(title = paste0(tag, " | ", term_name))

    ggsave(
      file.path(GSEA_DIR, "figs",
        paste0("enrichment_", tag, "_", i, ".png")),
      p,
      width = 8,
      height = 6
    )
  }
}

# ENRICHR

run_enrichr <- function(tag, genes) {

  if (!length(genes)) return(NULL)

  setEnrichrSite("FlyEnrichr")

  dbs <- grep("^GO_", listEnrichrDbs()$libraryName, value=TRUE)

  enrich <- enrichr(genes, dbs)

  for (d in names(enrich)) {
    fwrite(as.data.table(enrich[[d]]),
           file.path(GSEA_DIR,"tables",paste0(tag,"_",d,".csv")))
  }

  enrich
}

# NETWORK

run_network <- function(enrich_bp, tag) {

  if (is.null(enrich_bp) || nrow(enrich_bp)==0) return()

  top <- enrich_bp %>%
    arrange(Adjusted.P.value) %>%
    head(12)

  edge_list <- lapply(seq_len(nrow(top)), function(i) {

    genes <- unlist(strsplit(top$Genes[i], ";"))

    genes <- genes[genes != "" & !is.na(genes)]

    if (!length(genes)) return(NULL)

    data.frame(
      from = rep(top$Term[i], length(genes)),
      to   = genes,
      stringsAsFactors=FALSE
    )
  })

  edge_list <- edge_list[!sapply(edge_list,is.null)]
  if (!length(edge_list)) return()

  edges <- do.call(rbind, edge_list)

  g <- graph_from_data_frame(edges)

  p <- ggraph(g, layout="fr") +
    geom_edge_link(alpha=0.3) +
    geom_node_point(size=3) +
    geom_node_text(aes(label=name), repel=TRUE, size=3) +
    theme_void()

  ggsave(file.path(GSEA_DIR,"figs",paste0("network_",tag,".png")),
         p,width=10,height=8)
}

# MAIN

run_analysis <- function(tag) {

  message("\nRUNNING: ", tag)

  res <- load_de(tag)
  ranks <- build_ranks(res)

  pathways <- build_go_bp_sets(ranks)

  if (!is.null(pathways) && length(pathways)) {

    fg <- fgseaMultilevel(pathways, ranks)
    fg_dt <- as.data.table(fg)

    if (nrow(fg_dt) > 0) {

  fwrite(fg_dt,
    file.path(GSEA_DIR,"tables",paste0("fgsea_",tag,".csv"))
  )

  save_fgsea(fg_dt, tag)
  save_enrichment_plots(pathways, ranks, fg_dt, tag)
}
  }

  sig <- get_sig(res)
  sig_symbols <- fbgn_to_symbol(sig)

  enrich <- run_enrichr(tag, sig_symbols)

  if (!is.null(enrich[["GO_Biological_Process_2018"]])) {
    run_network(enrich[["GO_Biological_Process_2018"]], tag)
  }
}

# RUN ALL

for (t in CONTRASTS) {
  tryCatch(run_analysis(t),
    error=function(e) message("ERROR: ", e$message))
}

message("\nDONE → ", GSEA_DIR)