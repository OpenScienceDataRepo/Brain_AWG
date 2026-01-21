# Data Collection #
The following processed proteomic files from OSD-514 were used for analysis:
1. GLDS-514_proteomics_TMTc.tar.gz
2. GLDS-514_proteomics_TMTb.tar.gz
3. GLDS-514_proteomics_TMTa.tar.gz

# Processing Workflow #
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
install.packages("pheatmap")
install.packages("matrixStats")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(limma)
library(pheatmap)
library(matrixStats)
library(tibble)

setwd("/Users/mariannavarro/Desktop/ADBR Mito/OSD-514/Proteomic")
list.files()

# READ NASA METADATA
meta <- read_tsv(
  "a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt",
  show_col_types = FALSE
) %>%
  transmute(
    sample    = str_replace_all(`Sample Name`, " ", "_"),
    tmt_run   = factor(`Parameter Value[Run Number]`),
    label     = Label,
    condition = if_else(str_detect(sample, "^Earth"), "Earth", "Flight"),
    sex = case_when(
      str_detect(sample, "_M") ~ "Male",
      str_detect(sample, "_F") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(condition = factor(condition),
         sex = factor(sex))

tmt_runs <- levels(meta$tmt_run)

# FILE DISCOVERY
files_tbl <- map_dfr(tmt_runs, function(run) {
  tibble(
    tmt_run = run,
    file = list.files(run, pattern = "\\.txt$", full.names = TRUE)
  )
})

# parse sample+tmt_run from filename
files_tbl2 <- files_tbl %>%
  mutate(
    base   = basename(file) %>% str_remove("\\.txt$"),
    tmt_run = str_extract(base, "TMT[a-zA-Z0-9]+"),
    sample  = str_remove(base, "_TMT[a-zA-Z0-9]+$")
  )

sample_files <- files_tbl2 %>% filter(!str_detect(base, "^pool"))
pool_files   <- files_tbl2 %>% filter(str_detect(base, "^pool"))

# READ PROTEIN FILE FUNCTION
read_protein_file <- function(file) {

  lines <- readLines(file, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  if (length(lines) < 2) stop("File has no data rows: ", file)

  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]

  body_txt <- paste(lines[-1], collapse = "\n")
  df <- read.table(
    text = body_txt,
    sep = "\t",
    header = FALSE,
    quote = "\"",
    fill = TRUE,
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(df) > length(header)) {
    header <- c(header, paste0("X_extra_", seq_len(ncol(df) - length(header))))
  } else if (ncol(df) < length(header)) {
    header <- header[seq_len(ncol(df))]
  }
  names(df) <- header

  if (!("Accession" %in% names(df))) {
    stop("Expected 'Accession' column not found in: ", file,
         "\nColumns were:\n", paste(names(df), collapse = ", "))
  }

  abundance_cols <- names(df)[
    str_detect(names(df), regex("abundance", ignore_case = TRUE)) &
      !str_detect(names(df), regex("abundances\\s*count", ignore_case = TRUE))
  ]
  if (length(abundance_cols) == 0) {
    stop("No abundance column found in: ", file,
         "\nColumns were:\n", paste(names(df), collapse = ", "))
  }

  df %>%
    transmute(
      protein_id = as.character(Accession),
      intensity  = if (length(abundance_cols) == 1) {
        as.numeric(.data[[abundance_cols[1]]])
      } else {
        rowMeans(across(all_of(abundance_cols), as.numeric), na.rm = TRUE)
      }
    )
}

# READ SAMPLE DATA
sample_long <- sample_files %>%
  left_join(meta, by = c("sample" = "sample", "tmt_run" = "tmt_run")) %>%
  filter(!is.na(condition)) %>%
  mutate(data = map(file, read_protein_file)) %>%
  select(sample, tmt_run, condition, sex, data) %>%
  unnest(data) %>%
  group_by(protein_id, sample, tmt_run, condition, sex) %>%
  summarise(intensity = mean(intensity, na.rm = TRUE), .groups = "drop")

# BUILD EXPRESSION MATRIX
expr_mat <- sample_long %>%
  select(protein_id, sample, intensity) %>%
  pivot_wider(names_from = sample, values_from = intensity) %>%
  column_to_rownames("protein_id") %>%
  as.matrix()
mode(expr_mat) <- "numeric"

# READ POOL CHANNELS
pool_data <- pool_files %>%
  mutate(data = map(file, read_protein_file)) %>%
  unnest(data) %>%
  group_by(tmt_run, protein_id) %>%
  summarise(intensity = mean(intensity, na.rm = TRUE), .groups = "drop")

# INTERNAL REFERENCE SCALING (IRS)
expr_irs <- expr_mat

for (run in tmt_runs) {
  run_samples <- meta %>%
    filter(tmt_run == run) %>%
    pull(sample) %>%
    intersect(colnames(expr_irs))

  if (length(run_samples) == 0) next

  pool_vec <- pool_data %>%
    filter(tmt_run == run) %>%
    select(protein_id, intensity) %>%
    deframe()

  if (length(pool_vec) == 0) next

  common <- intersect(rownames(expr_irs), names(pool_vec))
  if (length(common) == 0) next

  expr_irs[common, run_samples] <- sweep(expr_irs[common, run_samples], 1, pool_vec[common], "/")
}

# LOG TRANSFORM & FILTER

expr_log <- log2(expr_irs)

# Remove proteins with >30% missing
keep_proteins <- rowSums(!is.na(expr_log)) >= 0.7 * ncol(expr_log)
expr_filt <- expr_log[keep_proteins, ]

# Remove proteins with zero variance
expr_filt <- expr_filt[rowVars(expr_filt, na.rm = TRUE) > 0, ]

# Remove samples with all NA
expr_filt <- expr_filt[, colSums(!is.na(expr_filt)) > 0]

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

# LIMMA DIFFERENTIAL ANALYSIS

common_samples <- intersect(colnames(expr_filt), meta$sample)
expr_limma <- expr_filt[, common_samples]

meta_limma <- meta %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

stopifnot(ncol(expr_limma) == nrow(meta_limma))

# Simple design: condition only
design <- model.matrix(~ condition, meta_limma)

fit <- lmFit(expr_limma, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = "conditionFlight", number = Inf, adjust.method = "BH")
results_df <- results %>% rownames_to_column("protein_id")
write_csv(results_df, "GLDS-514_TMT_limma_results.csv")

# VOLCANO PLOT

results_df <- results_df %>%
  mutate(sig = adj.P.Val < 0.05 & abs(logFC) > 1)

png("Volcano_Flight_vs_Earth.png", 800, 600)
ggplot(results_df, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = sig), alpha = 0.6) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal()
dev.off()

# HEATMAP (TOP 50)

top_proteins <- results_df %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(protein_id)

heat_mat <- expr_limma[top_proteins, ]
annotation_col <- meta_limma %>%
  select(sample, condition, sex, tmt_run) %>%
  column_to_rownames("sample")

png("Heatmap_Top50.png", 900, 900)
pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE
)
dev.off()
```
