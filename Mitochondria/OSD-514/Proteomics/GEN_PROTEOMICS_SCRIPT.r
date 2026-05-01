# LIMMA DIFFERENTIAL ANALYSIS

# PACKAGES
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("limma", ask = FALSE, update = FALSE)

install.packages(c("pheatmap", "matrixStats", "ggrepel"), dependencies = TRUE)

suppressPackageStartupMessages({
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(limma)
library(pheatmap)
library(matrixStats)
library(data.table)
library(fgsea)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(GO.db)
library(enrichR) 
library(tidyr) 
library(igraph) 
library(ggraph)
})

select <- dplyr::select
filter <- dplyr::filter

# WORKING DIRECTORY
setwd("/Volumes/Marians_SSD/ADBR_Mito/OSD-514/Proteomics")

# Output dirs
OUT_DIR <- "/Volumes/Marians_SSD/ADBR_Mito/OSD-514/Proteomics/RESULTS_OSD514"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figs"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(OUT_DIR, "figs")
TAB_DIR <- file.path(OUT_DIR, "tables")

# LOAD EXPRESSION MATRIX (FROM QC SCRIPT)
expr_mat <- read.csv("TMT_expression_matrix.csv", row.names = 1)
expr_mat <- as.matrix(expr_mat)
mode(expr_mat) <- "numeric"

# LOAD METADATA
meta <- read_tsv(
  "a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt",
  show_col_types = FALSE
) %>%
  transmute(
    sample = str_replace_all(`Sample Name`, " ", "_"),
    tmt_run = factor(`Parameter Value[Run Number]`),
    condition = case_when(
      str_detect(`Sample Name`, "^Earth") ~ "Earth",
      str_detect(`Sample Name`, "^SF1g") ~ "SF1g",
      str_detect(`Sample Name`, "^SFug") ~ "SFug",
      TRUE ~ NA_character_
    ),
    sex = case_when(
      str_detect(`Sample Name`, "_M") ~ "Male",
      str_detect(`Sample Name`, "_F") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition))

meta$condition <- factor(meta$condition, levels = c("Earth", "SF1g", "SFug"))

# ALIGN METADATA WITH EXPRESSION MATRIX
common_samples <- intersect(colnames(expr_mat), meta$sample)

expr_limma <- expr_mat[, common_samples]
meta_limma <- meta %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

stopifnot(all(colnames(expr_limma) == meta_limma$sample))

# DESIGN MATRIX
design <- model.matrix(~0 + condition, meta_limma)
colnames(design) <- levels(meta_limma$condition)

contrast_matrix <- makeContrasts(
  SF1g_vs_Earth = SF1g - Earth,
  SFug_vs_Earth = SFug - Earth,
  SF1g_vs_SFug  = SF1g - SFug,
  levels = design
)

# LIMMA fitting
fit <- lmFit(expr_limma, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Volcano Plots
save_volcano <- function(fit_obj, contrast_name,
                         pval_cut = 0.05,
                         logfc_cut = 1,
                         top_n = 5) {

  tt <- topTable(fit_obj, coef = contrast_name,
                 number = Inf, adjust.method = "BH") %>%
    rownames_to_column("protein_id") %>%
    mutate(
      sig = case_when(
        adj.P.Val < pval_cut & logFC >= logfc_cut  ~ "Up",
        adj.P.Val < pval_cut & logFC <= -logfc_cut ~ "Down",
        TRUE ~ "NotSig"
      )
    )

  top_labels <- bind_rows(
    tt %>% filter(sig == "Up") %>% arrange(desc(logFC)) %>% slice_head(n = top_n),
    tt %>% filter(sig == "Down") %>% arrange(logFC) %>% slice_head(n = top_n)
  )

  p <- ggplot(tt, aes(x = logFC, y = -log10(P.Value), color = sig)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-logfc_cut, logfc_cut),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_cut),
               linetype = "dashed") +
    geom_text_repel(
      data = top_labels,
      aes(label = protein_id),
      max.overlaps = 50,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      Down = "blue",
      Up = "red",
      NotSig = "grey70"
    )) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      x = "Log2 Fold Change",
      y = expression(-log[10](italic(p)))
    )

  ggsave(file.path(FIG_DIR, paste0("Volcano_", contrast_name, ".png")),
       p, width = 8, height = 6, dpi = 300)
}

# Save volcano plots
save_volcano(fit2, "SF1g_vs_Earth")
save_volcano(fit2, "SFug_vs_Earth")
save_volcano(fit2, "SF1g_vs_SFug")

results_SF1g_vs_Earth <- topTable(fit2, coef="SF1g_vs_Earth", number=Inf, adjust.method="BH") %>%
  rownames_to_column("protein_id")

results_SFug_vs_Earth <- topTable(fit2, coef="SFug_vs_Earth", number=Inf, adjust.method="BH") %>%
  rownames_to_column("protein_id")

results_SF1g_vs_SFug <- topTable(fit2, coef="SF1g_vs_SFug", number=Inf, adjust.method="BH") %>%
  rownames_to_column("protein_id")

top_proteins <- results_SF1g_vs_Earth %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(protein_id)

heat_mat <- expr_limma[top_proteins, ]

ann <- meta_limma %>%
  dplyr::select(sample, condition, sex, tmt_run) %>%   # <-- add dplyr::
  column_to_rownames("sample")

# Order columns
ord <- order(ann$condition, ann$sex, ann$tmt_run)
heat_mat <- heat_mat[, ord]
ann <- ann[ord, , drop=FALSE]

heat_z <- t(scale(t(heat_mat)))
heat_z[is.na(heat_z)] <- 0

png(file.path(FIG_DIR, "Heatmap_Top50.png"), 2200, 1500, res=200)
pheatmap(
  heat_z,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann,
  show_rownames = TRUE,
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA,
  main = "Top 50 Differential Proteins across Earth, SF1g, SFug"
)
dev.off()

# Save LIMMA Results Tables

write.csv(results_SF1g_vs_Earth,
          file.path(TAB_DIR, "Limma_SF1g_vs_Earth_results.csv"),
          row.names = FALSE)

write.csv(results_SFug_vs_Earth,
          file.path(TAB_DIR, "Limma_SFug_vs_Earth_results.csv"),
          row.names = FALSE)

write.csv(results_SF1g_vs_SFug,
          file.path(TAB_DIR, "Limma_SF1g_vs_SFug_results.csv"),
          row.names = FALSE)

# Save Only Sig Results Tables

sig_SF1g_vs_Earth <- results_SF1g_vs_Earth %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

write.csv(sig_SF1g_vs_Earth,
          file.path(TAB_DIR, "Limma_SF1g_vs_Earth_significant.csv"),
          row.names = FALSE)

sig_SFug_vs_Earth <- results_SFug_vs_Earth %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

write.csv(sig_SFug_vs_Earth,
          file.path(TAB_DIR, "Limma_SFug_vs_Earth_significant.csv"),
          row.names = FALSE)

sig_SF1g_vs_SFug <- results_SF1g_vs_SFug %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

write.csv(sig_SF1g_vs_SFug,
          file.path(TAB_DIR, "Limma_SF1g_vs_SFug_significant.csv"),
          row.names = FALSE)

# Quick Summary (In-Console)

cat("SF1g vs Earth:", nrow(sig_SF1g_vs_Earth), "significant proteins\n")
cat("SFug vs Earth:", nrow(sig_SFug_vs_Earth), "significant proteins\n")
cat("SF1g vs SFug:", nrow(sig_SF1g_vs_SFug), "significant proteins\n")

## GSEA

# PATHS

BASE_DIR <- "/Volumes/Marians_SSD/ADBR_Mito/OSD-514/Proteomics"

GSEA_DIR <- file.path(BASE_DIR, "GSEA")
TBL_DIR  <- file.path(BASE_DIR, "RESULTS_OSD514", "tables")

dir.create(file.path(GSEA_DIR, "tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(GSEA_DIR, "figs"), recursive=TRUE, showWarnings=FALSE)

# PARAMETERS

PADJ_CUT <- 0.05
LFC_CUT  <- 1
CONTRASTS <- c("SF1g_vs_Earth","SFug_vs_Earth","SF1g_vs_SFug")

# HELPERS

`%||%` <- function(a,b) if (!is.null(a) && length(a)) a else b
wrap_terms <- function(x,w=40) stringr::str_wrap(x,w)

# LOAD DE

load_de <- function(tag) {
  f <- file.path(TBL_DIR, paste0("Limma_", tag, "_results.csv"))
  if (!file.exists(f)) stop("Missing DE file: ", f)

  df <- read.csv(f, stringsAsFactors=FALSE, check.names=FALSE)

  rownames(df) <- df$protein_id
  df
}

# RANKS

build_ranks <- function(res) {

  ids <- rownames(res)

  pv  <- res$P.Value %||% res$adj.P.Val
  lfc <- res$logFC %||% 0

  pv[!is.finite(pv)] <- 1
  lfc[!is.finite(lfc)] <- 0

  score <- sign(lfc) * (-log10(pmax(pv, 1e-300)))

  keep <- is.finite(score) & nzchar(ids)
  ranks <- tapply(score[keep], ids[keep], mean)

  sort(ranks, decreasing = TRUE)
}

# GO BP SETS (RESTORED)

build_go_bp_sets <- function(ranks) {

  genes <- names(ranks)

  m <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys=genes,
    keytype="UNIPROT",
    columns=c("GOALL","ONTOLOGYALL")
  )

  m <- m[m$ONTOLOGYALL=="BP", ]

  if (nrow(m)==0) return(NULL)

  paths <- split(m$UNIPROT, m$GOALL)
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

# SIGNIFICANT PROTEINS

get_sig <- function(res) {
  res <- res[
    !is.na(res$adj.P.Val) &
    res$adj.P.Val < PADJ_CUT &
    !is.na(res$logFC) &
    abs(res$logFC) >= LFC_CUT,
  ]
  unique(rownames(res))
}

# SYMBOL CONVERSION

uniprot_to_symbol <- function(ids) {
  tbl <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = ids,
    keytype = "UNIPROT",
    columns = "SYMBOL"
  )

  tbl <- tbl[!duplicated(tbl$UNIPROT), ]
  map <- setNames(tbl$SYMBOL, tbl$UNIPROT)

  out <- map[ids]
  out[is.na(out)] <- ids[is.na(out)]

  out
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
  ranks_uniprot <- build_ranks(res)   # keep original UniProt IDs

  # Build GO sets using clean UniProt ranks
  pathways <- build_go_bp_sets(ranks_uniprot)

  if (is.null(pathways) || length(pathways) < 10) {
    message("Skipping ", tag, ": too few pathways")
    return(NULL)
  }

  # convert to symbols for fgsea
  ranks_named <- ranks_uniprot
  names(ranks_named) <- uniprot_to_symbol(names(ranks_uniprot))

  fg <- fgseaMultilevel(pathways, ranks_uniprot)  # use UniProt ranks here
  fg_dt <- as.data.table(fg)

  if (nrow(fg_dt) > 0) {
    fwrite(fg_dt, file.path(GSEA_DIR, "tables", paste0("fgsea_", tag, ".csv")))
    save_fgsea(fg_dt, tag)
    save_enrichment_plots(pathways, ranks_uniprot, fg_dt, tag)
  }

  sig <- get_sig(res)
  sig_symbols <- uniprot_to_symbol(sig)

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