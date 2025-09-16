# Procedure Overview

To analyze the 5xFAD data Study, differential expression was performed with DESeq2 in R across cortex and hippocampus at 4–18 months, followed by enrichment tests focused on mitochondrial and ROS pathways. Results showed a consistent oxidative-stress signal (NOX2–NRF2 axis) with only marginal shifts in OXPHOS genes at the pathway level. We therefore interpret the dataset as indicating clear ROS/stress activation in 5xFAD but limited evidence for robust respiratory-chain dysfunction at the transcript level.

# Metadata Preparation

The matrix is available on the GEO of the study but no metadata for R to read/analyze. Therefore, the csv had to be prepared before DESeq2 could be run. The below code was used:

```
# 0) Folder
setwd("~/Desktop/ADBR Mito/5xFAD/")  # <-- adjust to location of files
counts_file <- "counts.csv"
matrix_file <- "matrix.csv"
stopifnot(file.exists(counts_file), file.exists(matrix_file))

# 1) Read counts (TAB-delimited!)
counts  <- read.delim(counts_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
samples <- colnames(counts)
if (!length(samples)) stop("counts.csv has 0 sample columns. Is it TAB-delimited?")

# 2) Try normal GEO parse first
sm <- tryCatch(
  read.delim(matrix_file, header = FALSE, quote = "\"", check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) NULL
)

src_vec <- NULL
if (!is.null(sm)) {
  colnames(sm)[1] <- "TAG"
  idx_src <- which(sm$TAG == "!Sample_source_name_ch1")
  if (length(idx_src)) {
    src_row <- sm[idx_src, -1, drop = FALSE]
    src_vec <- trimws(gsub("^\"|\"$", "", as.character(unlist(src_row))))
  }
}

# 3) Fallback parser for files that use literal "\t" instead of real tabs
if (is.null(src_vec) || length(src_vec) <= 1) {
  lines <- readLines(matrix_file, warn = FALSE)
  line  <- lines[grep("^!Sample_source_name_ch1", lines, useBytes = TRUE)]
  if (!length(line)) stop("No !Sample_source_name_ch1 found in matrix.csv")
  # Remove the tag, normalize separators: convert literal \t to real tabs, strip quotes
  line2 <- sub("^!Sample_source_name_ch1", "", line, perl = TRUE)
  line2 <- gsub("^\\s+|^\\t+", "", line2, perl = TRUE)
  line2 <- gsub('\\"', '"', line2, perl = TRUE)           # normalize escaped quotes
  line2 <- gsub('\\\\t', '\t', line2, perl = TRUE)        # literal \t -> real tab
  line2 <- gsub('"', '', line2, perl = TRUE)              # drop quotes
  src_vec <- strsplit(line2, "\t", fixed = TRUE)[[1]]
  src_vec <- trimws(src_vec)
}

# 4) Sanity check: counts columns must match series-matrix samples
if (length(src_vec) != length(samples)) {
  cat("counts columns:", length(samples), " | matrix source_name entries:", length(src_vec), "\n")
  stop("Mismatch: counts columns != matrix source_name entries. Confirm both files are from the same series/run.")
}

# 5) Build metadata aligned by order
metadata <- data.frame(source_name = src_vec, stringsAsFactors = FALSE)
rownames(metadata) <- samples  # align to counts column order

# 6) Parse source_name -> strain / tissue / age / sex
split_tokens <- strsplit(metadata$source_name, "_", fixed = TRUE)

is_num <- function(x) !is.na(suppressWarnings(as.numeric(x)))
extract_fields <- function(tok) {
  tok <- tok[nzchar(tok)]
  # Drop trailing numeric ID if present (e.g., ..._Female_430)
  if (length(tok) >= 1 && is_num(tok[length(tok)])) tok <- tok[-length(tok)]
  n <- length(tok)
  sex    <- if (n >= 1) tolower(tok[n])       else NA_character_
  age    <- if (n >= 2) tolower(tok[n-1])     else NA_character_
  tissue <- if (n >= 3) tolower(tok[n-2])     else NA_character_
  strain <- if (n >= 4) paste(tok[1:(n-3)], collapse = "_") else NA_character_
  c(strain = strain, tissue = tissue, age = age, sex = sex)
}
parts_mat <- t(vapply(split_tokens, extract_fields, FUN.VALUE = c(strain="", tissue="", age="", sex="")))

clean <- function(x) {
  x <- gsub("[ ;]+", "_", x)
  x <- gsub("__+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}
parts_mat[, "strain"] <- clean(gsub(";", "", parts_mat[, "strain"]))  # remove semicolons like "5xFAD;BL6"
parts_mat[, "tissue"] <- clean(parts_mat[, "tissue"])
parts_mat[, "age"]    <- clean(gsub("\\bmonths?\\b", "mon", parts_mat[, "age"]))
parts_mat[, "sex"]    <- clean(parts_mat[, "sex"])

sex    <- parts_mat[, "sex"];    sex[is.na(sex) | sex==""] <- "unknown"
cg_raw <- paste(parts_mat[, "strain"], parts_mat[, "tissue"], parts_mat[, "age"], sep = "_")
cg_raw <- toupper(clean(cg_raw)); cg_raw[is.na(cg_raw) | cg_raw=="" | grepl("^NA_", cg_raw)] <- "UNKNOWN_GROUP"

meta_out <- data.frame(
  sex = sex,
  condition_group = cg_raw,
  row.names = rownames(metadata),
  stringsAsFactors = FALSE
)

# 7) Write metadata.csv
write.csv(meta_out, "metadata.csv")

cat("✅ metadata.csv written with", nrow(meta_out), "samples and columns: sex, condition_group\n")
cat("   Example (first 5):\n")
print(head(cbind(sample = rownames(meta_out), source_name = metadata$source_name, meta_out), 5))
```

# DESeq2 Analysis
```
# ===================== DESeq2 for multi-age 5xFAD vs BL6 (age numeric; coef-based shrink) =====================
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(biomaRt)
})

# 0) Working directory
setwd("~/Desktop/ADBR Mito/5xFAD/")  # <-- change if needed

# 1) Load data (counts is TAB-delimited!)
counts   <- read.delim("counts.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv("metadata.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(metadata) == ncol(counts),
          "condition_group" %in% colnames(metadata),
          "sex" %in% colnames(metadata))

# 2) Parse condition_group = STRAIN_TISSUE_AGE -> fields
tok <- strsplit(as.character(metadata$condition_group), "_", fixed = TRUE)
get_or_na <- function(x,i) if (length(x) >= i) x[[i]] else NA_character_
metadata$strain_raw <- toupper(vapply(tok, get_or_na, "", i = 1))
metadata$tissue     <- tolower(vapply(tok, get_or_na, "", i = 2))
metadata$age_label  <- toupper(vapply(tok, get_or_na, "", i = 3))
metadata$age_num    <- suppressWarnings(as.numeric(sub("MON","", metadata$age_label, ignore.case = TRUE)))
stopifnot(!any(is.na(metadata$age_num)))  # require numeric ages

# Normalize strain to {BL6, 5XFAD}; set factors
metadata$strain <- ifelse(grepl("5XFAD", metadata$strain_raw, ignore.case = TRUE), "5XFAD", "BL6")
metadata$strain <- factor(metadata$strain, levels = c("BL6","5XFAD"))  # BL6 reference
metadata$tissue <- factor(metadata$tissue)
metadata$sex    <- factor(tolower(metadata$sex), levels = c("female","male","unknown"))

# Use centered numeric age for overall model (helps GLM stability)
metadata$age      <- as.numeric(scale(metadata$age_num, center = TRUE, scale = FALSE))
metadata$age_label <- paste0(metadata$age_num, "mon")  # clean label for within-age tests

# 3) Align column order
counts <- counts[, rownames(metadata)]

# 4) Integer counts + strip Ensembl version suffix
rownames(counts) <- sub("\\..*$", "", rownames(counts))  # ENSMUSG...".4" -> ENSMUSG...
counts <- round(as.matrix(counts))

# 5) (Optional) Mouse gene annotation
try({
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genes <- getBM(
    attributes = c("ensembl_gene_id","external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = rownames(counts),
    mart       = ensembl
  )
  sym_map <- setNames(genes$external_gene_name, genes$ensembl_gene_id)
  gene_symbols <- sym_map[rownames(counts)]
  gene_symbols[is.na(gene_symbols) | gene_symbols==""] <- rownames(counts)[is.na(gene_symbols) | gene_symbols==""]
  rownames(counts) <- make.unique(gene_symbols)
}, silent = TRUE)

# 6) Helper: volcano plot
make_volcano <- function(res, title, lfc_cutoff = 1, p_cutoff = 0.05) {
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  res$significance <- "Not significant"
  res$significance[res$padj < p_cutoff & res$log2FoldChange >  lfc_cutoff] <- "Upregulated"
  res$significance[res$padj < p_cutoff & res$log2FoldChange < -lfc_cutoff] <- "Downregulated"
  res$negLogP <- -log10(res$padj)
  ggplot(res, aes(x = log2FoldChange, y = negLogP, color = significance)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = c("Upregulated"="red","Downregulated"="blue","Not significant"="grey")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    labs(title = title, x = "log2(Fold Change)", y = "-log10(Adjusted p-value)") +
    theme_minimal()
}

# -------------------------------
# A) OVERALL MODEL (age numeric)
# -------------------------------
dds_all <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata[, c("sex","tissue","age","strain")],
  design    = ~ sex + tissue + age + strain
)
dds_all <- dds_all[rowSums(counts(dds_all)) > 10, ]
dds_all <- DESeq(dds_all)

# Shrink overall strain effect
rn <- resultsNames(dds_all)
coef_strain <- rn[grepl("^strain_.*_vs_", rn)]
stopifnot(length(coef_strain) == 1)
shrink_type_overall <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
res_overall <- lfcShrink(dds_all, coef = coef_strain, type = shrink_type_overall)

write.csv(as.data.frame(res_overall), "DESeq2_overall_5xFAD_vs_BL6.csv")
p_overall <- make_volcano(res_overall, "Overall 5xFAD vs BL6 (adjusted for sex+tissue+age)")
ggsave("Volcano_overall_5xFAD_vs_BL6.png", p_overall, width = 7, height = 6, dpi = 300)

# ---------------------------------------
# B) WITHIN-AGE & WITHIN-TISSUE CONTRASTS
#    (5xFAD vs BL6 at each age per tissue)
# ---------------------------------------
# Build a group factor combining strain.tissue.age; keep sex as covariate (not stratified)
metadata$group <- interaction(metadata$strain, metadata$tissue, metadata$age_label, drop = TRUE, sep = ".")
metadata$group <- droplevels(metadata$group)

dds_grp <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata[, c("sex","group")],
  design    = ~ sex + group
)
dds_grp <- dds_grp[rowSums(counts(dds_grp)) > 10, ]
dds_grp <- DESeq(dds_grp)  # fit size factors & dispersions once

tissues <- levels(metadata$tissue)
ages    <- sort(unique(metadata$age_label))
have_apeglm <- requireNamespace("apeglm", quietly = TRUE)
have_ashr   <- requireNamespace("ashr",   quietly = TRUE)

mk <- function(s,t,a) paste(s, t, a, sep = ".")
for (t in tissues) {
  for (a in ages) {
    g5  <- mk("5XFAD", t, a)
    gbl <- mk("BL6",   t, a)
    if (!all(c(g5, gbl) %in% levels(metadata$group))) {
      message("Skipping ", t, " ", a, " (one group missing)")
      next
    }

    # Relevel reference to BL6.tissue.age, reuse dispersions; only redo Wald test
    dds_ref <- dds_grp
    dds_ref$group <- droplevels(relevel(dds_ref$group, ref = gbl))
    dds_ref <- nbinomWaldTest(dds_ref)  # quick; no refit of size factors/dispersion

    coef_name <- paste0("group_", g5, "_vs_", gbl)
    tag <- paste0("DESeq2_5xFAD_vs_BL6_", t, "_", a)

    if (have_apeglm && coef_name %in% resultsNames(dds_ref)) {
      res <- lfcShrink(dds_ref, coef = coef_name, type = "apeglm")
    } else if (have_ashr) {
      res <- lfcShrink(dds_ref, contrast = c("group", g5, gbl), type = "ashr")
    } else {
      res <- lfcShrink(dds_ref, contrast = c("group", g5, gbl), type = "normal")
    }

    write.csv(as.data.frame(res), paste0(tag, ".csv"))
    p <- make_volcano(res, paste("5xFAD vs BL6 —", t, " ", a))
    ggsave(paste0("Volcano_5xFAD_vs_BL6_", t, "_", a, ".png"), p, width = 7, height = 6, dpi = 300)
    message("Wrote: ", tag, ".csv and Volcano_5xFAD_vs_BL6_", t, "_", a, ".png")
  }
}

message("Done:
 - Overall: DESeq2_overall_5xFAD_vs_BL6.csv + Volcano_overall_5xFAD_vs_BL6.png
 - Per age/tissue: DESeq2_5xFAD_vs_BL6_<tissue>_<age>.csv + Volcano_*.png")
```

# Pre-GSEA File Creation & Preparation

A ranked gene list for GSEA was created using the below code

```
# 0) Set working directory
setwd("~/Desktop/ADBR Mito/5xFAD/")  # <- change to the folder with your DESeq2_*.csv files

# 1) Load packages
suppressPackageStartupMessages({ library(dplyr); library(tools) })

# 2) Robust ranker: prefer 'stat' -> signed -log10(pvalue) -> log2FoldChange
rank_from_results <- function(df) {
  stopifnot(nrow(df) > 0)
  cn <- tolower(colnames(df))

  has_stat <- "stat" %in% cn
  has_p    <- "pvalue" %in% cn
  has_lfc  <- "log2foldchange" %in% cn

  if (has_stat) {
    score <- df[[which(cn == "stat")]]
  } else if (has_p && has_lfc) {
    lfc <- df[[which(cn == "log2foldchange")]]
    p   <- df[[which(cn == "pvalue")]]
    score <- sign(lfc) * (-log10(p))
  } else if (has_lfc) {
    score <- df[[which(cn == "log2foldchange")]]
  } else {
    stop("Could not find 'stat', 'pvalue'+'log2FoldChange', or 'log2FoldChange' in this table.")
  }

  tibble(Gene = rownames(df), Score = as.numeric(score)) %>%
    filter(!is.na(Gene), Gene != "", is.finite(Score)) %>%
    # if the same symbol appears multiple times, keep the strongest signal
    mutate(absScore = abs(Score)) %>%
    arrange(desc(absScore)) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    select(Gene, Score) %>%
    arrange(desc(Score))
}

# 3) Function to write .rnk from a single CSV
generate_rnk <- function(file_in, file_out) {
  res <- read.csv(file_in, row.names = 1, check.names = FALSE)
  rnk <- rank_from_results(res)
  write.table(rnk, file = file_out, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message("Saved: ", file_out, "  (", nrow(rnk), " genes)")
}

# 4) Find all result tables and convert
csvs <- list.files(pattern = "^DESeq2_.*\\.csv$", full.names = TRUE)
stopifnot(length(csvs) > 0)

for (f in csvs) {
  base <- file_path_sans_ext(basename(f))       # e.g., "DESeq2_5xFAD_vs_BL6_cortex_8mon"
  out  <- paste0("GSEA_ranks_", base, ".rnk")   # e.g., "GSEA_ranks_DESeq2_5xFAD_vs_BL6_cortex_8mon.rnk"
  tryCatch(
    generate_rnk(f, out),
    error = function(e) message("Skipping ", f, " -> ", conditionMessage(e))
  )
}

message("Done. .rnk files are ready for GSEA/fgsea.")
```

# GSEA Analysis and Results

GSEA reinforces the DESeq2 findings: a clear ROS/oxidative-stress signal at mid-disease with later, cortex-specific mitochondrial effects. In cortex at 12 months, fatty-acid oxidation is enriched (NES≈1.72, FDR≈0.006), and by 18 months there is strong enrichment for OXPHOS/oxidative metabolism (e.g., Complex I–IV NES≈2.60, FDR≈0.010; oxidative metabolism NES≈2.54, FDR≈0.012) In contrast, hippocampus and the overall analysis show little to no mito-pathway enrichment at standard FDR thresholds, indicating that mitochondrial changes are age-progressive and region-specific.
