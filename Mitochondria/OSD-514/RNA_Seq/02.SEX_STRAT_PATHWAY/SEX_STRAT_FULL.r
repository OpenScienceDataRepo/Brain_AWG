## LIBRARIES
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(fgsea)
  library(clusterProfiler)
  library(org.Dm.eg.db)
  library(GO.db)
})

## OUTPUT DIRS
OUT_DIR <- "/Volumes/Marians_SSD/ADBR_Mito/OSD-514/RNA_Seq/RESULTS_Sex_Strat"
dir.create(file.path(OUT_DIR, "figs"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive=TRUE, showWarnings=FALSE)

GSEA_DIR <- file.path("/Volumes/Marians_SSD/ADBR_Mito/OSD-514/RNA_Seq/GSEA_Sex_Strat")
dir.create(file.path(GSEA_DIR, "figs"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(GSEA_DIR, "tables"), recursive=TRUE, showWarnings=FALSE)

## 1. IMPORT FILES + METADATA

files <- list.files(
  path = "/Volumes/Marians_SSD/ADBR_Mito/OSD-514/RNA_Seq/RawCounts",
  pattern = "\\.genes\\.results$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(files) == 0) {
  stop("No *.genes.results files found. Make sure you've downloaded the 24 '...genes.results' files from OSDR.")
}
message("Found ", length(files), " *.genes.results files.")

b <- basename(files)

m <- stringr::str_match(
  b,
  "rna-seq_(SFug|SF1g|Earth)_([MF])(\\d).*_(CRRA\\d+)\\-[^_]+_(HV\\w+)_L(\\d)"
)

stopifnot(!any(is.na(m)))

cond_map <- c(
  SFug  = "SPACEFLIGHT_MICROGRAVITY",
  SF1g  = "SPACEFLIGHT_1G",
  Earth = "EARTH"
)

sex_map <- c(M = "MALE", F = "FEMALE")

condition_group <- cond_map[m[,2]]
sex <- sex_map[m[,3]]
replicate <- as.integer(m[,4])

samp <- paste(condition_group, sex, replicate)
names(files) <- samp

meta <- data.frame(
  condition_group = factor(condition_group,
                           levels=c("EARTH","SPACEFLIGHT_1G","SPACEFLIGHT_MICROGRAVITY")),
  sex = factor(sex, levels=c("FEMALE","MALE")),
  replicate = replicate,
  row.names = samp
)

stopifnot(all(names(files) == rownames(meta)))

## 2. IMPORT + QC FILTERING

txi <- tximport(files, type="rsem", countsFromAbundance="no")

txi$length[txi$length <= 0 | is.na(txi$length)] <- 1

keep <- rowSums(txi$counts >= 10) >= ceiling(ncol(txi$counts) * 0.2)

txi$counts <- txi$counts[keep, ]
txi$abundance <- txi$abundance[keep, ]
txi$length <- txi$length[keep, ]

message("Kept genes after QC filter: ", nrow(txi$counts))

## 3. DESEQ2 MODEL

meta$condition_group <- relevel(meta$condition_group, "EARTH")

design_formula <- ~ condition_group * sex

dds <- DESeqDataSetFromTximport(txi, meta, design_formula)
dds <- DESeq(dds)

saveRDS(dds, file.path(OUT_DIR, "tables/dds.rds"))

## 4. QC (PCA + DISPERSION)

vsd <- vst(dds, blind=TRUE)

pca_df <- plotPCA(vsd, intgroup=c("condition_group","sex"), returnData=TRUE)
percentVar <- round(100 * attr(pca_df,"percentVar"))

p <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color=condition_group, shape=sex), size=3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_classic()

ggsave(file.path(OUT_DIR,"figs/pca.png"), p, width=7, height=6)

png(file.path(OUT_DIR,"figs/dispersion.png"))
plotDispEsts(dds)
dev.off()

## 5. DIFFERENTIAL EXPRESSION + store results objects

library(ggrepel)

toSym <- function(ids) {
  ids_clean <- gsub("\\..*$", "", ids)
  syms <- mapIds(org.Dm.eg.db,
                 keys      = ids_clean,
                 column    = "SYMBOL",
                 keytype   = "FLYBASE",
                 multiVals = "first")
  ifelse(is.na(syms), ids, syms)
}

contrasts <- list(
  MG_vs_E  = c("SPACEFLIGHT_MICROGRAVITY", "EARTH"),
  G1_vs_E  = c("SPACEFLIGHT_1G",           "EARTH"),
  MG_vs_1G = c("SPACEFLIGHT_MICROGRAVITY", "SPACEFLIGHT_1G")
)

res_list <- list()

for (nm in names(contrasts)) {
  a <- contrasts[[nm]][1]; b <- contrasts[[nm]][2]
  res     <- results(dds, contrast=c("condition_group", a, b))
  res_shr <- lfcShrink(dds, contrast=c("condition_group", a, b), type="ashr")
  res_list[[nm]] <- list(res=res, res_shr=res_shr)
  message(nm, " done — ", sum(as.data.frame(res)$padj < 0.05, na.rm=TRUE), " sig genes")
}

## 6. VISUALS — volcano

save_volcano <- function(res, res_shr, tag, thr_p=0.05, thr_fc=1) {
  
  df <- as.data.frame(res)
  df$gene       <- rownames(df)
  df$log2FC_shr <- as.data.frame(res_shr)[df$gene, "log2FoldChange"]
  df$gene_symbol <- toSym(df$gene)
  
  df$group <- "NS"
  df$group[!is.na(df$padj) & df$padj < thr_p & df$log2FoldChange >  thr_fc] <- "Up"
  df$group[!is.na(df$padj) & df$padj < thr_p & df$log2FoldChange < -thr_fc] <- "Down"
  
  lab_set <- head(df$gene[order(df$padj)], 15)
  
  out <- file.path(OUT_DIR, "figs", paste0("volcano_", tag, ".png"))
  png(out, width=1800, height=1500, res=200, bg="white")
  print(
    ggplot(df, aes(x=log2FC_shr, y=-log10(padj), color=group)) +
      geom_point(alpha=0.7, size=1.6, na.rm=TRUE) +
      geom_point(data=subset(df, gene %in% lab_set),
                 shape=21, stroke=0.7, size=2.8, color="black") +
      geom_vline(xintercept=c(-thr_fc, thr_fc), linetype="dashed") +
      geom_hline(yintercept=-log10(thr_p),       linetype="dashed") +
      geom_text_repel(
        data=subset(df, gene %in% lab_set),
        aes(label=gene_symbol),
        size=3, max.overlaps=40, box.padding=0.4, min.segment.length=0
      ) +
      scale_color_manual(values=c("Down"="#2C7BB6","NS"="grey60","Up"="#D7191C")) +
      labs(x="Shrunken log2 fold change", y="-log10(FDR)",
           title=paste0("Volcano: ", gsub("_"," ", tag)),
           subtitle="Dashed: |log2FC|=1 and FDR=0.05") +
      theme_classic() + theme(legend.position="top")
  )
  dev.off()
  message("Saved: ", out)
}

## Call once per contrast — pull from res_list
save_volcano(res_list[["MG_vs_E"]]$res,  res_list[["MG_vs_E"]]$res_shr,  "MG_vs_EARTH")
save_volcano(res_list[["G1_vs_E"]]$res,  res_list[["G1_vs_E"]]$res_shr,  "1G_vs_EARTH")
save_volcano(res_list[["MG_vs_1G"]]$res, res_list[["MG_vs_1G"]]$res_shr, "MG_vs_1G")


get_res <- function(a, b) {
  results(dds, contrast=c("condition_group", a, b))
}

## 7. GO ORA
for (nm in names(contrasts)) {
  res <- get_res(contrasts[[nm]][1], contrasts[[nm]][2])
  sig <- rownames(res)[which(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1)]
  if (length(sig) < 10) { message("Skipping GO ORA for ", nm, " — too few sig genes"); next }
  ego <- enrichGO(
    gene          = sig,
    OrgDb         = org.Dm.eg.db,
    keyType       = "FLYBASE",
    ont           = "BP",
    pAdjustMethod = "BH"
  )
  fwrite(as.data.table(ego),
         file.path(OUT_DIR,"tables",paste0("GO_",nm,".csv")))
  message("GO ORA saved: ", nm)
}

## mito heatmaps
mito_fbgn_curated <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys    = "GO:0005739",          # mitochondrion (CC)
  keytype = "GOALL",
  columns = "FLYBASE"
)$FLYBASE
mito_fbgn_curated <- unique(mito_fbgn_curated[!is.na(mito_fbgn_curated)])
message("Curated mito genes: ", length(mito_fbgn_curated))

## 8. HEATMAPS

save_mito_heatmap <- function(tag, res_unshr, groups, mito_fbgn, max_rows=50, cap_z=2.5) {
  if (!exists("vsd") || !exists("meta")) { message("No vsd/meta; skipping heatmap for ", tag); return(invisible(NULL)) }
  mat <- assay(vsd)

  cols_use <- rownames(meta)[meta$condition_group %in% groups]
  if (!length(cols_use)) { message("No columns match groups for ", tag); return(invisible(NULL)) }

  ann <- meta[cols_use, c("condition_group","sex","replicate"), drop=FALSE]

  cond_levels <- c("EARTH","SPACEFLIGHT_1G","SPACEFLIGHT_MICROGRAVITY")
  ord <- order(
    factor(ann$condition_group, levels=cond_levels[cond_levels %in% groups]),
    ann$sex,
    ann$replicate
  )
  ann <- ann[ord, , drop=FALSE]
  cols_use <- rownames(ann)

  cond_short <- c(
    EARTH                    = "Earth",
    SPACEFLIGHT_1G           = "SF1g",
    SPACEFLIGHT_MICROGRAVITY = "SFug"
  )
  sex_short <- c(FEMALE="F", MALE="M")

  lab_col <- paste0(
    cond_short[as.character(ann$condition_group)], "_",
    sex_short[as.character(ann$sex)],
    ann$replicate
  )

  rdf <- as.data.frame(res_unshr); rdf$gene_id <- rownames(rdf)
  genes_in <- intersect(mito_fbgn, rownames(mat))
  if (!length(genes_in)) { message("Curated mito set has no overlap for ", tag); return(invisible(NULL)) }
  genes_rankable <- intersect(genes_in, rownames(rdf))
  ord_idx <- order(rdf[genes_rankable, "padj"], -abs(rdf[genes_rankable, "log2FoldChange"]))
  genes_plot <- genes_rankable[ord_idx][seq_len(min(length(genes_rankable), max_rows))]

  z <- t(scale(t(mat[genes_plot, cols_use, drop=FALSE])))
  z[is.na(z)] <- 0
  z[z >  cap_z] <-  cap_z
  z[z < -cap_z] <- -cap_z
  rownames(z)   <- toSym(rownames(z))
  colnames(z)   <- lab_col
  rownames(ann) <- lab_col
  ann$replicate <- NULL

  out <- file.path(OUT_DIR,"figs", paste0("heatmap_mito_", tag, ".png"))
  png(out, width=2200, height=1500, res=200, bg="white")
  pheatmap::pheatmap(z,
    scale="none",
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    annotation_col = ann[, c("condition_group","sex"), drop=FALSE],
    show_rownames=TRUE, show_colnames=TRUE,
    fontsize_row=8, fontsize_col=9,
    border_color=NA,
    main=paste0("Curated mitochondrial genes (", gsub("_"," ", tag), "); row-Z cap ±", cap_z,
                " (", nrow(z), " genes)")
  )
  dev.off()
  message("Saved: ", out)
}

save_mito_heatmap("MG_vs_EARTH",
  res_unshr = res_list[["MG_vs_E"]]$res,
  groups    = c("SPACEFLIGHT_MICROGRAVITY","EARTH"),
  mito_fbgn = mito_fbgn_curated, max_rows=50, cap_z=2.5)

save_mito_heatmap("1G_vs_EARTH",
  res_unshr = res_list[["G1_vs_E"]]$res,
  groups    = c("SPACEFLIGHT_1G","EARTH"),
  mito_fbgn = mito_fbgn_curated, max_rows=50, cap_z=2.5)

save_mito_heatmap("MG_vs_1G",
  res_unshr = res_list[["MG_vs_1G"]]$res,
  groups    = c("SPACEFLIGHT_MICROGRAVITY","SPACEFLIGHT_1G"),
  mito_fbgn = mito_fbgn_curated, max_rows=50, cap_z=2.5)

message("Saved focused PNGs to: ", file.path(OUT_DIR, "figs"))

## 9. GSEA (uses res_list so no get_res needed here)

build_go_bp <- function(ranks) {
  genes <- names(ranks)
  m <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys    = genes,
    keytype = "FLYBASE",
    columns = c("GOALL","ONTOLOGYALL")
  )
  m <- m[m$ONTOLOGYALL=="BP", ]
  m <- m[!is.na(m$GOALL), ]
  if (nrow(m)==0) return(NULL)
  paths <- split(m$FLYBASE, m$GOALL)
  paths <- lapply(paths, unique)
  lens  <- lengths(paths)
  paths[lens >= 10 & lens <= 500]
}

get_go_terms <- function(go_ids) {
  if (!length(go_ids)) return(character(0))
  tbl <- AnnotationDbi::select(GO.db, keys=go_ids, keytype="GOID", columns="TERM")
  tbl <- unique(tbl)
  setNames(tbl$TERM, tbl$GOID)
}

plot_fgsea_bar <- function(fg_dt, tag) {
  if (nrow(fg_dt)==0) return()
  fg_dt    <- fg_dt[order(padj)][1:min(15,nrow(fg_dt))]
  term_map <- get_go_terms(fg_dt$pathway)
  fg_dt$term  <- term_map[fg_dt$pathway]
  fg_dt$term[is.na(fg_dt$term)] <- fg_dt$pathway
  fg_dt$label <- stringr::str_wrap(fg_dt$term, 40)
  p <- ggplot(fg_dt, aes(reorder(label,NES), NES, fill=-log10(padj))) +
    geom_col() + coord_flip() + theme_minimal() +
    labs(title=paste("GSEA:", tag))
  ggsave(file.path(GSEA_DIR,"figs",paste0("fgsea_",tag,".png")), p, width=10, height=6)
}

plot_enrichment <- function(pathways, ranks, fg_dt, tag) {
  if (nrow(fg_dt)==0) return()
  top      <- fg_dt[order(padj)][1:min(5,nrow(fg_dt))]
  term_map <- get_go_terms(top$pathway)
  for (i in seq_len(nrow(top))) {
    pw    <- top$pathway[i]
    label <- term_map[pw]
    if (is.na(label)) label <- pw
    p <- fgsea::plotEnrichment(pathways[[pw]], ranks) +
      ggtitle(paste(tag, "|", label))
    ggsave(file.path(GSEA_DIR,"figs",paste0("enrichment_",tag,"_",i,".png")), p, width=8, height=6)
  }
}

for (nm in names(contrasts)) {
  message("\n[GSEA] ", nm)
  res    <- as.data.frame(res_list[[nm]]$res)
  ranks  <- res$log2FoldChange
  names(ranks) <- rownames(res)
  ranks  <- sort(ranks[is.finite(ranks)], decreasing=TRUE)

  pathways <- build_go_bp(ranks)
  if (is.null(pathways)) { message("No pathways: ", nm); next }

  fg    <- fgseaMultilevel(pathways, ranks)
  fg_dt <- as.data.table(fg)

  fwrite(fg_dt, file.path(GSEA_DIR,"tables",paste0("fgsea_",nm,".csv")))
  plot_fgsea_bar(fg_dt, nm)
  plot_enrichment(pathways, ranks, fg_dt, nm)
  message("DONE GSEA: ", nm)
}
