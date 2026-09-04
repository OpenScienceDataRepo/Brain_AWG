# OSD-514 SEX-STRATIFIED RNASeq

library(DESeq2)
library(tximport)
library(dplyr)

BASE_DIR <- "/Volumes/Marians_SSD/ADBR_Mito/OSD-514/RNA_Seq"
RAW_DIR <- "/Volumes/Marians_SSD/ADBR_Mito/OSD-514/RNA_Seq/RawCounts"
META_FILE <- file.path(BASE_DIR, "Collapsed_Counts/metadata_from_filenames.csv")

OUT_DIR <- file.path(BASE_DIR, "RESULTS_OSD514")

# 1. LOAD METADATA
meta <- read.csv(META_FILE, row.names = 1)

meta$condition_group <- factor(meta$condition_group,
                               levels = c("EARTH",
                                          "SPACEFLIGHT_1G",
                                          "SPACEFLIGHT_MICROGRAVITY"))

meta$sex <- factor(meta$sex)

# 2. TXIMPORT
files <- list.files(RAW_DIR, pattern = "genes.results", full.names = TRUE)
names(files) <- gsub(".*/|\\.genes.results", "", files)

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

txi$length[txi$length == 0] <- 1

# align metadata
meta <- meta[colnames(txi$counts), ]

# 3. HELPER FUNCTION
run_deseq <- function(meta_sub, txi_sub) {

  dds <- DESeqDataSetFromTximport(
    txi_sub,
    colData = meta_sub,
    design = ~ condition_group
  )

  dds <- DESeq(dds)

  list(
    MG_vs_EARTH = results(dds, contrast = c("condition_group",
                                            "SPACEFLIGHT_MICROGRAVITY",
                                            "EARTH")),
    G1_vs_EARTH = results(dds, contrast = c("condition_group",
                                            "SPACEFLIGHT_1G",
                                            "EARTH")),
    MG_vs_G1 = results(dds, contrast = c("condition_group",
                                         "SPACEFLIGHT_MICROGRAVITY",
                                         "SPACEFLIGHT_1G"))
  )
}

extract_summary <- function(res) {
  df <- as.data.frame(res)
  df <- df[!is.na(df$padj), ]
  sig <- df[df$padj < 0.05, ]

  data.frame(
    sig_genes = nrow(sig),
    upregulated = sum(sig$log2FoldChange > 0),
    downregulated = sum(sig$log2FoldChange < 0)
  )
}

# 4. SPLIT DATA (THIS IS THE KEY FIX)
meta_male <- meta[meta$sex == "MALE", ]
meta_female <- meta[meta$sex == "FEMALE", ]

txi_male <- txi
txi_male$counts <- txi$counts[, rownames(meta_male)]
txi_male$length <- txi$length[, rownames(meta_male)]

txi_female <- txi
txi_female$counts <- txi$counts[, rownames(meta_female)]
txi_female$length <- txi$length[, rownames(meta_female)]

# 5. RUN ANALYSIS PER SEX
male_res <- run_deseq(meta_male, txi_male)
female_res <- run_deseq(meta_female, txi_female)

# 6. BUILD FIGURE 3A TABLE
summary_table <- data.frame(
  comparison = c(
    "MG_vs_EARTH_MALE",
    "1G_vs_EARTH_MALE",
    "MG_vs_G1_MALE",
    "MG_vs_EARTH_FEMALE",
    "1G_vs_EARTH_FEMALE",
    "MG_vs_G1_FEMALE"
  ),

  rbind(
    extract_summary(male_res$MG_vs_EARTH),
    extract_summary(male_res$G1_vs_EARTH),
    extract_summary(male_res$MG_vs_G1),
    extract_summary(female_res$MG_vs_EARTH),
    extract_summary(female_res$G1_vs_EARTH),
    extract_summary(female_res$MG_vs_G1)
  )
)

# 7. SAVE OUTPUT
write.csv(summary_table,
          file.path(OUT_DIR, "FIGURE3A_SEX_STRATIFIED_SUMMARY.csv"),
          row.names = FALSE)

cat("\nDONE ✔ Figure 3A-style results saved\n")
