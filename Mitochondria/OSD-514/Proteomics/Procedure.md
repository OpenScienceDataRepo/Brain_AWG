# Data Collection #
The following processed proteomic files from OSD-514 were used for analysis:
1. GLDS-514_proteomics_TMTc.tar.gz
2. GLDS-514_proteomics_TMTb.tar.gz
3. GLDS-514_proteomics_TMTa.tar.gz

# Processing Workflow #
```
# INSTALL & LOAD PACKAGES
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
install.packages(c("pheatmap","matrixStats","ggrepel"))

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
setwd("/Users/mariannavarro/Desktop/ADBR Mito/OSD-514/Proteomic")

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

# FILE DISCOVERY
files_tbl <- map_dfr(tmt_runs, function(run) {
  tibble(tmt_run = run, file = list.files(run, pattern = "\\.txt$", full.names = TRUE))
})

files_tbl2 <- files_tbl %>%
  mutate(
    base    = basename(file) %>% str_remove("\\.txt$"),
    tmt_run = str_extract(base, "TMT[a-zA-Z0-9]+"),
    sample  = str_remove(base, "_TMT[a-zA-Z0-9]+$")
  )

sample_files <- files_tbl2 %>% filter(!str_detect(base, "^pool"))
pool_files   <- files_tbl2 %>% filter(str_detect(base, "^pool"))

# FUNCTION TO READ PROTEIN FILES
read_protein_file <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  df <- read.table(text = paste(lines[-1], collapse="\n"),
                   sep="\t", header=FALSE, quote="\"", fill=TRUE,
                   stringsAsFactors = FALSE, check.names=FALSE)
  if(ncol(df) != length(header)) header <- c(header, paste0("X_extra_", seq_len(ncol(df)-length(header))))
  names(df) <- header
  abundance_cols <- names(df)[str_detect(names(df), regex("abundance", ignore_case = TRUE)) &
                                !str_detect(names(df), regex("abundances\\s*count", ignore_case = TRUE))]
  df %>% transmute(protein_id = as.character(Accession),
                   intensity = if(length(abundance_cols)==1) as.numeric(.data[[abundance_cols[1]]])
                   else rowMeans(across(all_of(abundance_cols)), na.rm=TRUE))
}

# READ SAMPLE DATA
sample_long <- sample_files %>%
  left_join(meta, by = c("sample","tmt_run")) %>%
  filter(!is.na(condition)) %>%
  mutate(data = map(file, read_protein_file)) %>%
  select(sample, tmt_run, condition, sex, label, data) %>%
  unnest(data) %>%
  group_by(protein_id, sample, tmt_run, condition, sex, label) %>%
  summarise(intensity = mean(intensity, na.rm=TRUE), .groups="drop")

# BUILD EXPRESSION MATRIX
expr_mat <- sample_long %>%
  select(protein_id, sample, intensity) %>%
  pivot_wider(names_from=sample, values_from=intensity) %>%
  column_to_rownames("protein_id") %>%
  as.matrix()
mode(expr_mat) <- "numeric"

# READ POOL DATA AND IRS NORMALIZATION
pool_data <- pool_files %>% mutate(data = map(file, read_protein_file)) %>% unnest(data) %>%
  group_by(tmt_run, protein_id) %>% summarise(intensity=mean(intensity, na.rm=TRUE), .groups="drop")

expr_irs <- expr_mat
for(run in tmt_runs){
  run_samples <- intersect(meta %>% filter(tmt_run==run) %>% pull(sample), colnames(expr_irs))
  pool_vec <- pool_data %>% filter(tmt_run==run) %>% select(protein_id,intensity) %>% deframe()
  common <- intersect(rownames(expr_irs), names(pool_vec))
  if(length(run_samples)>0 & length(common)>0){
    expr_irs[common, run_samples] <- sweep(expr_irs[common, run_samples], 1, pool_vec[common], "/")
  }
}

# LOG TRANSFORM & FILTER
expr_log <- log2(expr_irs)
keep_proteins <- rowSums(!is.na(expr_log)) >= 0.7*ncol(expr_log)
expr_filt <- expr_log[keep_proteins, ]
expr_filt <- expr_filt[rowVars(expr_filt, na.rm=TRUE)>0, ]
expr_filt <- expr_filt[, colSums(!is.na(expr_filt))>0]

# PCA
pca_mat <- expr_filt[rowVars(expr_filt, na.rm = TRUE) > 1e-6, ]
pca_mat[is.na(pca_mat)] <- 0  # replace remaining NAs with 0 for PCA

pca <- prcomp(t(pca_mat), center = TRUE, scale = FALSE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(meta %>% select(sample, condition, tmt_run), by = "sample")

png("PCA_TMT_samples.png", 800, 600)
ggplot(pca_df, aes(PC1, PC2, color = condition, shape = tmt_run)) +
  geom_point(size = 4) +
  theme_minimal()
dev.off()

# --- LIMMA DIFFERENTIAL ANALYSIS WITH FLIGHT SUBTYPES ---

# Reorder meta_limma to match columns of expr_limma
meta_limma <- meta %>%
  mutate(
    condition2 = case_when(
      str_detect(sample, "^Earth") ~ "Earth",
      str_detect(sample, "^SF1g")  ~ "SF1g",
      str_detect(sample, "^SFug")  ~ "SFug"
    ),
    condition2 = factor(condition2, levels = c("Earth","SF1g","SFug"))
  )

# Match sample names exactly
common_samples <- intersect(colnames(expr_filt), meta_limma$sample)
expr_limma <- expr_filt[, common_samples]
meta_limma <- meta_limma %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

stopifnot(all(colnames(expr_limma) == meta_limma$sample))  # should be TRUE

# Verify alignment
all(colnames(expr_limma) == meta_limma$sample)  # should return TRUE
stopifnot(all(colnames(expr_limma) == meta_limma$sample))

design <- model.matrix(~0 + condition2, meta_limma)
colnames(design) <- levels(meta_limma$condition2)

contrast_matrix <- makeContrasts(
  SF1g_vs_Earth = SF1g - Earth,
  SFug_vs_Earth = SFug - Earth,
  SF1g_vs_SFug  = SF1g - SFug,
  levels = design
)
meta_limma <- meta %>%
  mutate(
    condition2 = case_when(
      str_detect(sample, "^Earth") ~ "Earth",
      str_detect(sample, "^SF1g")  ~ "SF1g",
      str_detect(sample, "^SFug")  ~ "SFug"
    ),
    condition2 = factor(condition2, levels = c("Earth","SF1g","SFug"))
  )
table(meta_limma$condition2)

# LIMMA fitting
fit <- lmFit(expr_limma, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Volcano Plots
save_volcano <- function(fit_obj, contrast_name, pval_cut=0.05, logfc_cut=1, top_n=5){
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  
  # Get top table
  tt <- topTable(fit_obj, coef=contrast_name, number=Inf, adjust.method="BH") %>%
    rownames_to_column("protein_id")
  
  # Add significance and direction
  tt <- tt %>%
    mutate(
      sig = case_when(
        adj.P.Val < pval_cut & logFC >= logfc_cut  ~ "Up",
        adj.P.Val < pval_cut & logFC <= -logfc_cut ~ "Down",
        TRUE ~ "NotSig"
      )
    )
  
  # Label top N upregulated and top N downregulated proteins
  top_up <- tt %>% filter(sig=="Up") %>% arrange(desc(logFC)) %>% slice_head(n=top_n)
  top_down <- tt %>% filter(sig=="Down") %>% arrange(logFC) %>% slice_head(n=top_n)
  top_labels <- bind_rows(top_up, top_down)
  
  # Volcano plot
  p <- ggplot(tt, aes(x=logFC, y=-log10(adj.P.Val), color=sig, label=protein_id)) +
    geom_point(alpha=0.7) +
    scale_color_manual(values=c("Down"="blue", "Up"="red", "NotSig"="grey")) +
    geom_text_repel(data = top_labels, max.overlaps=50) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value")
  
  # Save plot
  ggsave(paste0("Volcano_", contrast_name, ".png"), p, width=8, height=6, dpi=200)
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
  select(sample, condition2, sex, tmt_run) %>%
  column_to_rownames("sample")

# Order columns
ord <- order(ann$condition2, ann$sex, ann$tmt_run)
heat_mat <- heat_mat[, ord]
ann <- ann[ord, , drop=FALSE]

heat_z <- t(scale(t(heat_mat)))
heat_z[is.na(heat_z)] <- 0

png("Heatmap_Top50_Earth_SF1g_SFug.png", 2200, 1500, res=200)
pheatmap(
  heat_z,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann,
  show_rownames = FALSE,
  fontsize_col = 9,
  border_color = NA,
  main = "Top 50 Differential Proteins across Earth, SF1g, SFug"
)
dev.off()

```
