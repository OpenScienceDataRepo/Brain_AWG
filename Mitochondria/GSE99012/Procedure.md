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

```
# Setup Working Directory
setwd("~/Desktop/ADBR Mito/GSE99012")  # change if needed

# Load Packages
pkgs_cran <- c("ggplot2","ggrepel","pheatmap","dplyr","readxl")
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
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  library(GO.db)
})

# --- Output Directory ---
OUT_DIR <- "RESULTS_GSE99012"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# --- Load DEG Table (Table S1) ---
deg_table <- read_excel("DEG Sekiya et al. (2018).xlsx", sheet = "Table S1")
deg_table <- deg_table[-c(1,2), ]  # Remove metadata/header rows

# --- Rename columns ---
colnames(deg_table)[1:8] <- c("GeneSymbol","Symbol","human_ortholog","log2FC","t","PValue","FDR","Comparison")

# --- Convert numeric columns ---
deg_table$log2FC <- as.numeric(deg_table$log2FC)
deg_table$FDR    <- as.numeric(deg_table$FDR)
deg_table$GeneSymbol <- as.character(deg_table$GeneSymbol)
deg_table$Symbol <- as.character(deg_table$Symbol)
deg_table$Comparison <- as.character(deg_table$Comparison)

# --- Curated Drosophila Mitochondrial Genes (FBgn IDs via GO) ---
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

# --- Volcano Plot Function ---
save_volcano_mito <- function(df, tag, mito_genes, thr_FDR=0.05, thr_LFC=1) {
  
  # Skip if no numeric values
  if(all(is.na(df$log2FC)) | all(is.na(df$FDR))) {
    message("Skipping ", tag, ": no numeric log2FC/FDR")
    return(NULL)
  }
  
  # Assign significance groups
  df$group <- "NS"
  df$group[df$FDR < thr_FDR & df$log2FC >=  thr_LFC] <- "Up (sig)"
  df$group[df$FDR < thr_FDR & df$log2FC <= -thr_LFC] <- "Down (sig)"
  df$group <- factor(df$group, levels=c("Down (sig)","NS","Up (sig)"))
  
  # Label only significant mitochondrial genes
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

# --- Key Comparisons: Only conditions vs Control ---
comparisons_to_plot <- c(
  "eGRL_TREM2-WT/TYROBP vs eGRL_Control",
  "eGRL_TREM2-R47H/TYROBP vs eGRL_Control",
  "eGRL_Aß42 vs eGRL_Control",
  "eGRL_Aß42/TREM2-WT/TYROBP vs eGRL_Control",
  "eGRL_Aß42/TREM2-R47H/TYROBP vs eGRL_Control"
)

# --- Generate Volcano Plots for Each ---
for (comp in comparisons_to_plot) {
  df_sub <- deg_table[deg_table$Comparison == comp, ]
  save_volcano_mito(df_sub, tag=comp, mito_genes=mito_genes)
}

message("✅ Mitochondrial-focused volcano plots saved to: ", OUT_DIR)
```
