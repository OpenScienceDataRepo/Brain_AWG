# INSTALL & LOAD PACKAGES
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("limma", ask = FALSE, update = FALSE)

install.packages(c("pheatmap","matrixStats","ggrepel"), dependencies = TRUE)

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(limma)
library(pheatmap)
library(matrixStats)
library(tibble)

# SET WORKING DIRECTORY
setwd("/Volumes/Marians_SSD/ADBR_Mito/OSD-514/Proteomics")

# READ METADATA
meta <- read_tsv("a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt",
                 show_col_types = FALSE) %>%
  transmute(
    sample    = str_replace_all(`Sample Name`, " ", "_"),
    tmt_run   = factor(`Parameter Value[Run Number]`),
    label     = Label,
    condition = if_else(str_detect(sample, "^Earth"), "Earth", "Flight"),
    sex       = case_when(
      str_detect(sample, "_M") ~ "Male",
      str_detect(sample, "_F") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(condition = factor(condition),
         sex       = factor(sex))

tmt_runs <- levels(meta$tmt_run)

# BUILD EXPRESSION MATRIX
expr_mat <- read.csv("TMT_expression_matrix.csv", row.names = 1)
expr_mat <- as.matrix(expr_mat)
mode(expr_mat) <- "numeric"

# LOG TRANSFORM & FILTER
expr_log <- log2(expr_mat)

keep_proteins <- rowSums(!is.na(expr_log)) >= 0.7 * ncol(expr_log)
expr_filt <- expr_log[keep_proteins, ]
expr_filt <- expr_filt[rowVars(expr_filt, na.rm=TRUE) > 0, ]
expr_filt <- expr_filt[, colSums(!is.na(expr_filt)) > 0]

pca_mat <- expr_filt[rowVars(expr_filt, na.rm = TRUE) > 1e-6, ]
pca_mat[is.na(pca_mat)] <- 0

pca <- prcomp(t(pca_mat), center = TRUE, scale = FALSE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(meta %>% select(sample, condition, tmt_run), by = "sample")

png("PCA_TMT_samples.png", 800, 600)
ggplot(pca_df, aes(PC1, PC2, color = condition, shape = tmt_run)) +
  geom_point(size = 4) +
  theme_minimal()
dev.off()