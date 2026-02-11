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
# Set working directory
setwd("/Volumes/Marians SSD/ADBR Mito/PRJNA747152")

# Load packages
library(dplyr)
library(readxl)
library(tidyr)
library(pheatmap)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(GO.db)

# Create output directory
OUT_DIR <- "RESULTS_PRJNA747152"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Detect all DEG .xlsx files
deg_files <- list.files(pattern="^DEG.*\\.xlsx$", full.names=TRUE)
names(deg_files) <- gsub(".xlsx","",basename(deg_files))

# Load DEG tables (only FlyBase ID / gene symbol and log2FC)
deg_list <- lapply(deg_files, function(f){
  df <- read_excel(f)
  df <- as.data.frame(df)
  colnames(df)[1:2] <- c("GeneSymbol","log2FC")
  df$GeneSymbol <- as.character(df$GeneSymbol)
  df$log2FC <- as.numeric(df$log2FC)
  df$Comparison <- gsub(".xlsx","",basename(f))
  df[, c("GeneSymbol","log2FC","Comparison")]
})
deg_table <- bind_rows(deg_list)

# Define mitochondrial GO terms
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

go2genes <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = all_go,
  keytype = "GO",
  columns = "FLYBASE"
)
mito_genes <- unique(na.omit(go2genes$FLYBASE))
message("Number of mitochondrial genes detected: ", length(mito_genes))

# Map gene symbols to FlyBase IDs
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

# Heatmap function
save_heatmap_mito <- function(deg_table, mito_genes, tag, topN=40, cap_z=2.5){
  mat_fc <- deg_table %>%
    dplyr::filter(FBgn %in% mito_genes) %>%
    tidyr::pivot_wider(names_from=Comparison, values_from=log2FC, values_fill=0)

  if(nrow(mat_fc) == 0){
    message("No mito genes found for heatmap: ", tag)
    return(NULL)
  }

  # Use only numeric columns for max_abs_LFC
  numeric_cols <- sapply(mat_fc, is.numeric)
  mat_fc$max_abs_LFC <- apply(abs(mat_fc[, numeric_cols, drop=FALSE]), 1, max, na.rm=TRUE)
  mat_fc_top <- mat_fc %>% arrange(desc(max_abs_LFC)) %>% slice_head(n=min(topN, nrow(mat_fc)))

  row_labels <- mat_fc_top$GeneSymbol
  plot_mat <- as.matrix(mat_fc_top[, numeric_cols & colnames(mat_fc_top) != "max_abs_LFC", drop=FALSE])

  # Row-wise Z-score normalization
  z <- t(scale(t(plot_mat)))
  z[is.na(z)] <- 0
  z[z > cap_z] <- cap_z
  z[z < -cap_z] <- -cap_z

  # Save heatmap
  png(file.path(OUT_DIR, paste0("heatmap_", tag, ".png")), width=1800, height=1500, res=200)
  pheatmap(z,
           labels_row=row_labels,
           show_rownames=TRUE,
           cluster_rows=TRUE,
           cluster_cols=TRUE,
           main=paste0("Top ", nrow(z), " mito genes (row-Z cap ±", cap_z, ")"),
           fontsize_row=7, fontsize_col=9,
           border_color=NA)
  dev.off()
  message("Heatmap saved: ", tag)
}

# Generate heatmap
save_heatmap_mito(deg_table_mapped, mito_genes, tag="log2FC_only_mito")
```
