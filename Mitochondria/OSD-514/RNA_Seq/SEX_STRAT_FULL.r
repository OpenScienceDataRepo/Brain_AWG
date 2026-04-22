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

files <- list.files(".", pattern="\\.genes\\.results$", full.names=TRUE, recursive=TRUE)
stopifnot(length(files) > 0)

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

## 5. DIFFERENTIAL EXPRESSION

get_res <- function(a,b) {
  results(dds, contrast=c("condition_group", a, b))
}

contrasts <- list(
  MG_vs_E  = c("SPACEFLIGHT_MICROGRAVITY","EARTH"),
  G1_vs_E  = c("SPACEFLIGHT_1G","EARTH"),
  MG_vs_1G = c("SPACEFLIGHT_MICROGRAVITY","SPACEFLIGHT_1G")
)

for (nm in names(contrasts)) {

  res <- get_res(contrasts[[nm]][1], contrasts[[nm]][2])
  df <- as.data.frame(res)
  df$gene <- rownames(df)

  write.csv(df,
            file.path(OUT_DIR,"tables",paste0("DE_",nm,".csv")),
            row.names=FALSE)

  message(nm, " sig genes: ", sum(df$padj < 0.05, na.rm=TRUE))
}

## 6. VISUALS (VOLCANO + HEATMAP)

for (nm in names(contrasts)) {

  res <- get_res(contrasts[[nm]][1], contrasts[[nm]][2])
  df <- as.data.frame(res)
  df$gene <- rownames(df)

  df$group <- "NS"
  df$group[df$padj < 0.05 & df$log2FoldChange > 1] <- "Up"
  df$group[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"

  p <- ggplot(df, aes(log2FoldChange, -log10(padj), color=group)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept=c(-1,1), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    theme_classic() +
    labs(title=nm)

  ggsave(file.path(OUT_DIR,"figs",paste0("volcano_",nm,".png")), p)

  top <- rownames(df[order(df$padj),])[1:40]
  mat <- assay(vsd)[top, ]
  mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0

  ann <- meta[colnames(mat), c("condition_group","sex")]

  png(file.path(OUT_DIR,"figs",paste0("heatmap_",nm,".png")), width=1000, height=800)
  pheatmap(mat, annotation_col=ann, show_rownames=FALSE, main=nm)
  dev.off()
}

## 7. GO ENRICHMENT (ORA)

for (nm in names(contrasts)) {

  res <- get_res(contrasts[[nm]][1], contrasts[[nm]][2])

  sig <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]

  if (length(sig) < 10) next

  ego <- enrichGO(
    gene=sig,
    OrgDb=org.Dm.eg.db,
    keyType="FLYBASE",
    ont="BP",
    pAdjustMethod="BH"
  )

  fwrite(as.data.table(ego),
         file.path(OUT_DIR,"tables",paste0("GO_",nm,".csv")))
}

## 8. GSEA

make_ranks <- function(res) {
  r <- res$log2FoldChange
  names(r) <- rownames(res)
  r <- r[is.finite(r)]
  sort(r, decreasing=TRUE)
}

get_go_terms <- function(go_ids) {
  if (!length(go_ids)) return(character(0))

  tbl <- AnnotationDbi::select(
    GO.db,
    keys=go_ids,
    keytype="GOID",
    columns="TERM"
  )

  tbl <- unique(tbl)
  setNames(tbl$TERM, tbl$GOID)
}

build_go_bp <- function(ranks) {

  genes <- names(ranks)

  m <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys=genes,
    keytype="FLYBASE",
    columns=c("GOALL","ONTOLOGYALL")
  )

  m <- m[m$ONTOLOGYALL=="BP", ]
  m <- m[!is.na(m$GOALL), ]

  if (nrow(m)==0) return(NULL)

  paths <- split(m$FLYBASE, m$GOALL)
  paths <- lapply(paths, unique)

  lens <- lengths(paths)
  paths[lens >= 10 & lens <= 500]
}

plot_fgsea_bar <- function(fg_dt, tag) {

  if (nrow(fg_dt)==0) return()

  fg_dt <- fg_dt[order(padj)][1:min(15,nrow(fg_dt))]

  term_map <- get_go_terms(fg_dt$pathway)

  fg_dt$term <- term_map[fg_dt$pathway]
  fg_dt$term[is.na(fg_dt$term)] <- fg_dt$pathway

  fg_dt$label <- stringr::str_wrap(fg_dt$term, 40)

  p <- ggplot(fg_dt, aes(reorder(label,NES), NES, fill=-log10(padj))) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(title=paste("GSEA:",tag))

  ggsave(file.path(GSEA_DIR,"figs",paste0("fgsea_",tag,".png")),
         p,width=10,height=6)
}

plot_enrichment <- function(pathways,ranks,fg_dt,tag){

  if (nrow(fg_dt)==0) return()

  top <- fg_dt[order(padj)][1:min(5,nrow(fg_dt))]

  term_map <- get_go_terms(top$pathway)

  for(i in seq_len(nrow(top))) {

    pw <- top$pathway[i]
    label <- term_map[pw]
    if (is.na(label)) label <- pw

    p <- fgsea::plotEnrichment(pathways[[pw]], ranks) +
      ggtitle(paste(tag,"|",label))

    ggsave(file.path(GSEA_DIR,"figs",
                     paste0("enrichment_",tag,"_",i,".png")),
           p,width=8,height=6)
  }
}

for (nm in names(contrasts)) {

  message("\n[GSEA] ", nm)

  res <- get_res(contrasts[[nm]][1], contrasts[[nm]][2])

  ranks <- make_ranks(res)

  pathways <- build_go_bp(ranks)

  if (is.null(pathways)) {
    message("No pathways: ", nm)
    next
  }

  fg <- fgseaMultilevel(pathways, ranks)
  fg_dt <- as.data.table(fg)

  fwrite(fg_dt,
         file.path(GSEA_DIR,"tables",paste0("fgsea_",nm,".csv")))

  plot_fgsea_bar(fg_dt, nm)
  plot_enrichment(pathways, ranks, fg_dt, nm)

  message("DONE GSEA: ", nm)
}

cat("\nPIPELINE COMPLETE → ", OUT_DIR, "\n")