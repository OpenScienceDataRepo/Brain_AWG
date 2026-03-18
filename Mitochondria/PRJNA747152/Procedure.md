# FASTQ -> Count Matrices

## 1. FASTQ -> FASTQC

The FASTQ files for this dataset were downloaded from [this directory](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA747152). They were then converted to Count Matrices using the [FastQC package](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

To install the FastQC Package, I used the following code:
```
# changing directory to where my installation is
cd /Volumes/Marians\ SSD

# unzipping and moving the installation
unzip fastqc_v0.12.1.zip
mv FastQC ~/FastQC
cd ~/FastQC

# make fastqc executable
chmod 755 fastqc

# testing fastqc (if it launches, then it's ready to go)
./fastqc

# make FastQC available globally (using symlink)
sudo ln -s /Users/mariannavarro/FastQC/fastqc /usr/local/bin/fastqc

# test
fastqc --help
```

Once that was set up, I ran fastqc on the PRJNA747152 FastQ files using the following code:
```
fastqc /Volumes/Marians\ SSD/ADBR\ Mito/PRJNA747152/FASTQ/* -o /Volumes/Marians\ SSD/ADBR\ Mito/PRJNA747152
```

## 2. FASTQC -> BAM
I moved the resulting FastQC files to a folder called "FastQC". I then used the [HISAT2 package](https://daehwankimlab.github.io/hisat2/download/) to align the reads to the reference genome.

Note that some packages have to be installed and initialized if you don't have it already:
  [homebrew](https://brew.sh/)
  [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)
  [hisat2](https://daehwankimlab.github.io/hisat2/download/)
  [samtools](https://www.htslib.org/)

The following bash file was created for simplicity, titled "run_alignment_safe.sh":
```
#!/bin/bash

# Run this script with: bash run_alignment_safe.sh

# Paths
FASTQ_DIR="/Volumes/Marians_SSD/ADBR_Mito/PRJNA747152/FASTQ"
OUTPUT_DIR="/Volumes/Marians_SSD/ADBR_Mito/PRJNA747152/BAM"
INDEX="/Volumes/Marians_SSD/dm6/genome"

# Activate conda environment
echo "Activating rnaseq_env..."
# Ensure conda.sh is sourced so 'conda activate' works
source /opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh
conda activate rnaseq_env

# Confirm binaries
echo "Using HISAT2 at: $(which hisat2)"
echo "Using SAMtools at: $(which samtools)"

# Count total FASTQs
TOTAL=$(ls "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | wc -l)
BAM_COUNT=0

# Loop through FASTQ files
for fq in "$FASTQ_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    bam="$OUTPUT_DIR/${base}.bam"

    if [ -f "$bam" ]; then
        echo "[$base] BAM exists, skipping."
        BAM_COUNT=$((BAM_COUNT + 1))
        continue
    fi

    echo "[$base] Aligning FASTQ..."
    hisat2 -x "$INDEX" -U "$fq" -S "$OUTPUT_DIR/${base}.sam"

    echo "[$base] Sorting SAM to BAM..."
    samtools sort -o "$bam" "$OUTPUT_DIR/${base}.sam"
    rm "$OUTPUT_DIR/${base}.sam"

    BAM_COUNT=$((BAM_COUNT + 1))
    echo "=== Progress: $BAM_COUNT / $TOTAL BAMs completed ==="
done

echo "All files processed."
date
```
## 3. Run featureCounts (Quantification)
Using the [ensembl genome annotation](https://www.ensembl.org/) file and the [featureCounts package](https://academic.oup.com/bioinformatics/article/30/7/923/232889), we converted the BAM files into counts.
```
conda activate rnaseq_env    # reactivate the environment
conda install -c bioconda subread    # install featurecounts if you don't have it already
which featureCounts    # double check featurecounts is installed

featureCounts \
  -s 2 \
  -a "/Volumes/Marians_SSD/ADBR_Mito/Drosophila_melanogaster.BDGP6.54.115.gtf" \
  -o "/Volumes/Marians_SSD/ADBR_Mito/PRJNA747152/counts.txt" \
  /Volumes/Marians_SSD/ADBR_Mito/PRJNA747152/BAM/*.bam

# clean and trim counts list
awk 'NR==2{for(i=7;i<=NF;i++){gsub(".*/","",$i); gsub(".bam","",$i)}} {printf "%s", $1; for(i=2;i<=6;i++) printf "\t%s", $i; for(i=7;i<=NF;i++) printf "\t%s", $i; printf "\n"}' counts.txt > clean_counts.txt
tail -n +2 clean_counts.txt | cut -f1,6,7- > counts_gene_length_clean.txt

# check result
head -n 5 counts_gene_length_clean.txt
```
# Mito-Specific DEG and Visualization
```
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

# 1. Paths & Read Files

counts_file <- "counts_gene_length_clean.txt"
meta_file   <- "meta.csv"

counts <- read.delim(counts_file, row.names=1, check.names=FALSE)
counts_numeric <- counts[, 7:ncol(counts)]  # numeric columns

meta_full <- read.csv(meta_file, stringsAsFactors=FALSE)
sample_names <- colnames(counts_numeric)

# Subset metadata to match counts
meta <- meta_full[meta_full$Run %in% sample_names, ]
rownames(meta) <- meta$Run

# 2. Define condition groups

meta$condition_group <- sapply(meta$treatment, function(x) {
  if (x == "0 Gy gamma rays")        return("Control")
  if (x == "0.4 Gy gamma rays")      return("LowGamma")
  if (x == "10 Gy gamma rays")       return("HighGamma")
  if (x == "0.4+10Gy gamma rays")    return("ComboGamma")
  return(NA)
})
meta$condition_group <- factor(meta$condition_group,
                               levels=c("Control","LowGamma","HighGamma","ComboGamma"))

# Optional: sex as factor
if ("sex" %in% colnames(meta) && !is.factor(meta$sex)) meta$sex <- factor(meta$sex)

# Optional: batch (flowcell + lane)
if (!("batch" %in% colnames(meta)) && all(c("flowcell","lane") %in% colnames(meta))) {
  meta$batch <- factor(paste(meta$flowcell, meta$lane, sep="_"))
}

# 3. Prepare DESeq2 dataset

use_batch <- "batch" %in% colnames(meta) && nlevels(droplevels(meta$batch)) > 1
use_sex   <- "sex"   %in% colnames(meta) && nlevels(droplevels(meta$sex))   > 1

design_formula <- if (use_batch && use_sex) {
  ~ batch + sex + condition_group
} else if (use_batch) {
  ~ batch + condition_group
} else if (use_sex) {
  ~ sex + condition_group
} else {
  ~ condition_group
}
message("Design: ", deparse(design_formula))

dds <- DESeqDataSetFromMatrix(
  countData = counts_numeric,
  colData   = meta,
  design    = design_formula
)

# Optional: remove genes with zero counts
dds <- dds[rowSums(counts(dds)) > 0, ]

# Run DESeq
dds <- DESeq(dds)

# 4. DE Contrast Function

run_contrast_simple <- function(var, levelA, levelB, tag) {
  res <- results(dds, contrast=c(var, levelA, levelB), alpha=0.05)
  
  # Save unshrunk DE table
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  write.csv(res_df, file.path(OUT_DIR, "tables", paste0("DE_", tag, "_unshrunk.csv")), row.names=FALSE)
  
  # Shrink log2 fold change if coefficient exists
  coef_name <- paste0(var, "_", levelA, "_vs_", levelB)
  res_shr_df <- NULL
  if (coef_name %in% resultsNames(dds)) {
    res_shr <- lfcShrink(dds, coef=coef_name, type="apeglm")
    res_shr_df <- as.data.frame(res_shr)
    res_shr_df$gene <- rownames(res_shr)
    write.csv(res_shr_df, file.path(OUT_DIR, "tables", paste0("DE_", tag, "_LFCshrunk_apeglm.csv")), row.names=FALSE)
  }
  
  # DE summary
  sig <- subset(res_df, !is.na(padj) & padj < 0.05)
  sig_lfc1 <- subset(sig, abs(log2FoldChange) >= 1)
  message(sprintf("[%s] padj<0.05: %d; padj<0.05 & |LFC|>=1: %d", tag, nrow(sig), nrow(sig_lfc1)))
  
  invisible(list(res=res_df, res_shr=res_shr_df))
}

# 5. Run radiation contrasts

out_Low_vs_Control    <- run_contrast_simple("condition_group", "LowGamma",    "Control",      "Low_vs_Control")
out_High_vs_Control   <- run_contrast_simple("condition_group", "HighGamma",   "Control",      "High_vs_Control")
out_Combo_vs_Control  <- run_contrast_simple("condition_group", "ComboGamma", "Control",      "Combo_vs_Control")
out_Combo_vs_High     <- run_contrast_simple("condition_group", "ComboGamma", "HighGamma",    "Combo_vs_High")
out_Combo_vs_Low      <- run_contrast_simple("condition_group", "ComboGamma", "LowGamma",     "Combo_vs_Low")

# 6. Curated mitochondrial gene set

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

go2genes <- AnnotationDbi::select(org.Dm.eg.db, keys=all_go, keytype="GO", columns="FLYBASE")
mito_fbgn_curated <- unique(na.omit(go2genes$FLYBASE))
message("Curated mito set size (FlyBase IDs): ", length(mito_fbgn_curated))

# 7. Volcano plotting function (Option 1)

toSym <- function(x) {
  map_tbl <- AnnotationDbi::select(org.Dm.eg.db, keys=x, keytype="FLYBASE", columns="SYMBOL")
  sym <- setNames(map_tbl$SYMBOL, map_tbl$FLYBASE)
  y <- unname(sym[x]); y[is.na(y) | y==""] <- x[is.na(y) | y==""]
  y
}

save_volcano_mito <- function(res_unshr, res_shr, tag, mito_fbgn, thr_p=0.05, thr_fc=1) {
  vdf <- as.data.frame(res_unshr)
  vdf$gene_id <- rownames(vdf)
  
  # Use shrunk LFC if available, otherwise unshrunk
  if (!is.null(res_shr)) {
    vdf$log2FC_shr <- res_shr[rownames(vdf), "log2FoldChange"]
  } else {
    vdf$log2FC_shr <- vdf$log2FoldChange
  }
  
  vdf$gene_symbol <- toSym(vdf$gene_id)
  
  vdf$group <- "NS"
  vdf$group[vdf$padj < thr_p & vdf$log2FC_shr >= thr_fc] <- "Up (sig)"
  vdf$group[vdf$padj < thr_p & vdf$log2FC_shr <= -thr_fc] <- "Down (sig)"
  vdf$group <- factor(vdf$group, levels=c("Down (sig)","NS","Up (sig)"))
  
  lab_set <- intersect(mito_fbgn, vdf$gene_id)
  if (!length(lab_set)) lab_set <- head(vdf[order(vdf$padj), "gene_id"], 12)
  
  out <- file.path(OUT_DIR,"figs", paste0("volcano_mito_", tag, ".png"))
  png(out, width=1800, height=1500, res=200, bg="white")
  print(
    ggplot(vdf, aes(x=log2FC_shr, y=-log10(padj), color=group)) +
      geom_point(alpha=0.7, size=1.6, na.rm=TRUE) +
      geom_point(data=subset(vdf, gene_id %in% lab_set), shape=21, stroke=0.7, size=2.8, color="black") +
      geom_vline(xintercept=c(-thr_fc, thr_fc), linetype="dashed") +
      geom_hline(yintercept=-log10(thr_p), linetype="dashed") +
      ggrepel::geom_text_repel(
        data=subset(vdf, gene_id %in% lab_set),
        aes(label=gene_symbol),
        size=3, max.overlaps=40, box.padding=0.4, min.segment.length=0
      ) +
      scale_color_manual(values=c("Down (sig)"="#2C7BB6","NS"="grey60","Up (sig)"="#D7191C")) +
      labs(x="Log2 fold change", y="-log10(FDR)",
           title=paste0("Volcano (mitochondrial focus): ", gsub("_"," ", tag)),
           subtitle="Dashed: |log2FC|=1 and FDR=0.05. Labels = curated mito genes (GO).") +
      theme_classic() + theme(legend.position="top")
  )
  dev.off()
}

# 8. Generate volcano plots

save_volcano_mito(out_Low_vs_Control$res,    out_Low_vs_Control$res_shr,    "Low_vs_Control",    mito_fbgn_curated)
save_volcano_mito(out_High_vs_Control$res,   out_High_vs_Control$res_shr,   "High_vs_Control",   mito_fbgn_curated)
save_volcano_mito(out_Combo_vs_Control$res,  out_Combo_vs_Control$res_shr,  "Combo_vs_Control",  mito_fbgn_curated)
save_volcano_mito(out_Combo_vs_High$res,     out_Combo_vs_High$res_shr,     "Combo_vs_High",     mito_fbgn_curated)
save_volcano_mito(out_Combo_vs_Low$res,      out_Combo_vs_Low$res_shr,      "Combo_vs_Low",      mito_fbgn_curated)
```

# Citations
Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4

Homebrew. (2026). The Homebrew package manager for macOS / Linux. Retrieved from https://brew.sh/

Anaconda, Inc. (2026). Conda – Package, dependency, and environment management for any language. Retrieved from https://docs.conda.io/

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., … & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008
