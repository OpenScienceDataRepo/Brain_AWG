# Load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(tibble)

# Paths
BASE_DIR <- "/Volumes/Marians_SSD/ADBR_Mito/OSD-514/Proteomics"
TMT_DIR  <- file.path(BASE_DIR, "TMT_all")
meta_file <- file.path(BASE_DIR, "a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt")

# Load and clean metadata
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  transmute(
    raw_name = `Sample Name`,
    sample = str_replace_all(raw_name, " ", "_"),
    batch = `Parameter Value[Run Number]`,
    label = Label,
    condition = case_when(
      str_detect(sample, "^Earth") ~ "Earth",
      str_detect(sample, "^SF1g")  ~ "SF1g",
      str_detect(sample, "^SFug")  ~ "SFug",
      TRUE ~ NA_character_
    ),
    sex = case_when(
      str_detect(sample, "_M") ~ "Male",
      str_detect(sample, "_F") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition), !is.na(sex))

cat("Metadata samples:", nrow(meta), "\n")

# List TMT files
files <- list.files(TMT_DIR, pattern="\\.txt$", full.names = TRUE)

cat("Total files:", length(files), "\n")

# Function to read ONE sample file
read_tmt <- function(f){
  
  df <- read.delim(f, check.names = FALSE)
  
  # extract CLEAN sample name (REMOVE batch suffix HERE)
  sample_name <- basename(f) %>%
    str_remove("\\.txt$") %>%
    str_remove("_TMT[a-c]$")
  
  # remove pool files
  if(str_detect(sample_name, "^pool")) return(NULL)
  
  # find abundance column (should be exactly one)
  ab_col <- grep("^Abundance:", colnames(df), value = TRUE)
  
  if(length(ab_col) != 1){
    stop(paste("Unexpected abundance columns in:", sample_name))
  }
  
  tibble(
    protein = df$Accession,
    sample  = sample_name,
    value   = as.numeric(df[[ab_col]])
  ) %>%
    group_by(protein, sample) %>%
    summarise(value = median(value, na.rm = TRUE), .groups = "drop")
}

# Read all files
expr_long <- map_dfr(files, read_tmt)

cat("Samples detected:", length(unique(expr_long$sample)), "\n")

# Convert to matrix
expr_mat <- expr_long %>%
  pivot_wider(names_from = sample, values_from = value) %>%
  column_to_rownames("protein") %>%
  as.matrix()

# replace NA and log transform
expr_mat[is.na(expr_mat)] <- 0
expr_mat <- log2(expr_mat + 1)

cat("Expression matrix:", nrow(expr_mat), "proteins x", ncol(expr_mat), "samples\n")

# Align metadata with matrix
common <- intersect(meta$sample, colnames(expr_mat))

meta <- meta %>% filter(sample %in% common)
expr_mat <- expr_mat[, meta$sample, drop = FALSE]

cat("Samples after alignment:", length(common), "\n")

# Final sanity checks
cat("\nMissing in matrix:\n")
print(setdiff(meta$sample, colnames(expr_mat)))

cat("\nMissing in metadata:\n")
print(setdiff(colnames(expr_mat), meta$sample))

# Create output directory
out_dir <- file.path(BASE_DIR)

# Convert matrix to dataframe for writing
expr_df <- as.data.frame(expr_mat) %>%
  rownames_to_column("protein")

# Save file
out_file <- file.path(out_dir, "TMT_expression_matrix.csv")
write.csv(expr_df, out_file, row.names = FALSE)

cat("Saved combined matrix to:", out_file, "\n")