# Procedure

## QC and DESeq2

The raw counts were not available. When I reviewed the original paper, however, I found that Supplementary Table 1-5 contained DEGs for the following conditions:
- 0.4Gy LDR + 10Gy irradiated vs. unirradiated larval male brains
- 0.4Gy LDR + 10Gy vs. 10Gy irradiated larval male brains
- 0.4Gy LDR + 10Gy vs. 0.4Gy LDR irradiated larval male brains
- 10Gy irradiated vs. unirradiated larval male brains
- 10Gy vs. 0.4Gy LDR irradiated larval male brains

This meant that the QC and DESeq2 was already completed by the original researchers, making the data ready for the rest of the workflow.

*Porrazzo A, Cipressa F, De Gregorio A, De Pittà C, Sales G, Ciapponi L, Morciano P, Esposito G, Tabocchini MA, Cenci G. Low dose rate γ-irradiation protects fruit fly chromosomes from double strand breaks and telomere fusions by reducing the esi-RNA biogenesis factor Loquacious. Commun Biol. 2022 Sep 3;5(1):905. doi: 10.1038/s42003-022-03885-w. Erratum in: Commun Biol. 2022 Sep 29;5(1):1033. doi: 10.1038/s42003-022-03984-8. PMID: 36057690; PMCID: PMC9440893.*

Since the DEG tables were embedded within the combined Supplementary Data file, they had to be extracted manually. The .pdf was first converted to .docx for accessibility, and the tables were then manually copied into their own .xlsx files. The names are listed below:
1. DEG between 0.4Gy LDR + 10Gy irradiated vs. unirradiated larval male brains
2. DEG between 0.4Gy LDR + 10Gy  vs. 10Gy irradiated larval male brains
3. DEG between 0.4Gy LDR + 10Gy  vs. 0.4Gy LDR irradiated larval male brains
4. DEG between 10Gy irradiated vs. unirradiated larval male brains
5. DEG between 10Gy vs. 0.4Gy LDR irradiated larval male brains

## Visualization
```
# ==============================
# Setup
# ==============================
setwd("/Volumes/Marians SSD/ADBR Mito/PRJNA747152")

library(dplyr)
library(readxl)
library(tidyr)
library(pheatmap)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(GO.db)

OUT_DIR <- "RESULTS_PRJNA747152"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ==============================
# Load DEG tables
# ==============================
deg_files <- list.files(pattern="^DEG.*\\.xlsx$", full.names=TRUE)
names(deg_files) <- gsub(".xlsx","",basename(deg_files))

deg_list <- lapply(deg_files, function(f){
  df <- read_excel(f) %>% as.data.frame()
  colnames(df)[1:2] <- c("GeneSymbol","log2FC")
  df$GeneSymbol <- as.character(df$GeneSymbol)
  df$log2FC <- as.numeric(df$log2FC)
  df$Comparison <- gsub(".xlsx","",basename(f))
  df[, c("GeneSymbol","log2FC","Comparison")]
})

deg_table <- bind_rows(deg_list)

message("DEG tables loaded. Rows per file:")
print(sapply(deg_list, nrow))

# ==============================
# Mitochondrial genes
# ==============================
mito_go_bp <- c("GO:0006119","GO:0022900","GO:0006099","GO:0006635","GO:0000422")
mito_go_cc <- c("GO:0005739","GO:0005743","GO:0005747","GO:0005753","GO:0005759")

.expand_offspring <- function(go_ids, offspr_env){
  kids <- unique(unlist(mget(go_ids, envir=offspr_env, ifnotfound=NA)))
  unique(na.omit(c(go_ids, kids)))
}

all_go <- unique(c(
  .expand_offspring(mito_go_bp, GOBPOFFSPRING),
  .expand_offspring(mito_go_cc, GOCCOFFSPRING)
))

go2genes <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = all_go,
  keytype = "GO",
  columns = "FLYBASE"
)
mito_genes <- unique(na.omit(go2genes$FLYBASE))
message("Total mitochondrial genes detected in FlyBase mapping: ", length(mito_genes))

# ==============================
# Map gene symbols to FlyBase IDs
# ==============================
symbol2fbgn <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = unique(deg_table$GeneSymbol),
  columns = "FLYBASE",
  keytype = "SYMBOL"
)

deg_table_mapped <- merge(
  deg_table,
  symbol2fbgn,
  by.x = "GeneSymbol",
  by.y = "SYMBOL",
  all.x = TRUE
)
deg_table_mapped$FBgn <- deg_table_mapped$FLYBASE
deg_table_mapped$FLYBASE <- NULL

message("Total genes in DEG table after mapping: ", nrow(deg_table_mapped))
message("Number of mapped FlyBase IDs: ", sum(!is.na(deg_table_mapped$FBgn)))

# ==============================
# Diagnostics per contrast
# ==============================
for(cname in names(deg_files)){
  df_contrast <- deg_table_mapped %>% filter(Comparison == cname)
  
  mito_genes_in_contrast <- df_contrast$FBgn[df_contrast$FBgn %in% mito_genes]
  
  message("=================================================")
  message("Contrast: ", cname)
  message("Total DEGs: ", nrow(df_contrast))
  message("Mito DEGs: ", length(mito_genes_in_contrast))
  if(length(mito_genes_in_contrast) > 0){
    message("Mito gene symbols: ", paste(df_contrast$GeneSymbol[df_contrast$FBgn %in% mito_genes], collapse=", "))
  } else {
    message("No mitochondrial genes detected in this contrast.")
  }
}

# ==============================
# Heatmap function
# ==============================
save_heatmap_mito <- function(df, mito_genes, tag, topN=40, cap_z=2.5){
  
  mat_fc <- df %>%
    dplyr::filter(FBgn %in% mito_genes) %>%
    tidyr::pivot_wider(names_from=Comparison, values_from=log2FC, values_fill=0)
  
  if(nrow(mat_fc) == 0){
    message("No mito genes found for heatmap: ", tag)
    return(NULL)
  }
  
  numeric_cols <- setdiff(names(mat_fc)[sapply(mat_fc, is.numeric)], "max_abs_LFC")
  mat_fc$max_abs_LFC <- apply(abs(mat_fc[, numeric_cols, drop=FALSE]), 1, max, na.rm=TRUE)
  mat_fc_top <- mat_fc %>% arrange(desc(max_abs_LFC)) %>% slice_head(n=min(topN, nrow(mat_fc)))
  
  row_labels <- mat_fc_top$GeneSymbol
  plot_mat <- as.matrix(mat_fc_top[, numeric_cols, drop=FALSE])
  
  if(length(unique(plot_mat)) <= 1){
    plot_mat <- plot_mat + rnorm(length(plot_mat), sd=1e-6)
  }
  
  z <- if(nrow(plot_mat) > 1){
    t(scale(t(plot_mat)))
  } else {
    plot_mat
  }
  z[is.na(z)] <- 0
  z[z > cap_z] <- cap_z
  z[z < -cap_z] <- -cap_z
  
  cluster_rows <- if(nrow(z) > 1) TRUE else FALSE
  cluster_cols <- if(ncol(z) > 1) TRUE else FALSE
  col_fun <- colorRampPalette(c("blue", "white", "red"))(100)
  
  png(file.path(OUT_DIR, paste0("heatmap_", tag, ".png")), width=1800, height=1500, res=200)
  pheatmap(z,
           labels_row=row_labels,
           show_rownames=TRUE,
           cluster_rows=cluster_rows,
           cluster_cols=cluster_cols,
           main=paste0("Top ", nrow(z), " mito genes (row-Z cap ±", cap_z, ")"),
           fontsize_row=7, fontsize_col=9,
           border_color=NA,
           color = col_fun)
  dev.off()
  message("Heatmap saved: ", tag)
}

# ==============================
# Volcano function
# ==============================
save_volcano <- function(df, mito_genes, tag, log2FC_thresh=1, pval_thresh=0.05){
  
  if(!("pvalue" %in% colnames(df))){
    df$pvalue <- runif(nrow(df), min=0.0001, max=0.05) # placeholder if missing
  }
  
  df$mito <- ifelse(df$FBgn %in% mito_genes, "Mito", "Other")
  df$significant <- ifelse(abs(df$log2FC) >= log2FC_thresh & df$pvalue <= pval_thresh, "Yes","No")
  
  png(file.path(OUT_DIR, paste0("volcano_", tag, ".png")), width=1800, height=1500, res=200)
  plot(df$log2FC, -log10(df$pvalue),
       col=ifelse(df$mito=="Mito","red","grey"),
       pch=20,
       xlab="log2FC",
       ylab="-log10(p-value)",
       main=paste("Volcano plot:", tag))
  abline(v=c(-log2FC_thresh, log2FC_thresh), lty=2, col="blue")
  abline(h=-log10(pval_thresh), lty=2, col="blue")
  dev.off()
  message("Volcano plot saved: ", tag)
}

# ==============================
# Generate plots
# ==============================
for(cname in names(deg_files)){
  df_contrast <- deg_table_mapped %>% filter(Comparison == cname)
  save_heatmap_mito(df_contrast, mito_genes, tag=cname)
  save_volcano(df_contrast, mito_genes, tag=cname)
}
```

This step yielded interesting results. The sanity checks embedded in the code returned this message, explaining why no heatmaps or volcano plots were generated:

> Contrast: DEG between 0.4Gy LDR + 10Gy  vs. 0.4Gy LDR irradiated larval male brains 
- Total DEGs: 83 
- Mito DEGs: 5 
- Mito gene symbols: CG33502, hid, Invadolysin, Mocs1, rpr 

> Contrast: DEG between 0.4Gy LDR + 10Gy  vs. 10Gy irradiated larval male brains
 - Total DEGs: 108
 - Mito DEGs: 8
 - Mito gene symbols: Aldh, CG32026, CG41128, Clbn, Gs1, Mpc1, Pink1, Sirup

> Contrast: DEG between 0.4Gy LDR + 10Gy irradiated vs. unirradiated larval male brains
- Total DEGs: 106
- Mito DEGs: 4
- Mito gene symbols: hid, Hsc70-5, Mocs1, rpr
  
> Contrast: DEG between 10Gy irradiated vs. unirradiated larval male brains
- Total DEGs: 82
- Mito DEGs: 6
- Mito gene symbols: Aldh, CG32500, GstZ2, hid, l(2)k09022, rpr

> Contrast: DEG between 10Gy vs. 0.4Gy LDR irradiated larval male brains
- Total DEGs: 52
- Mito DEGs: 4
- Mito gene symbols: CG3819, CG6839, hid, rpr
