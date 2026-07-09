# LIMMA + GSEA + ENRICHR DIFFERENTIAL ANALYSIS (MS3FIX)

# PACKAGES
cran_pkgs <- c("dplyr", "tibble", "readr", "stringr", "ggplot2", "ggrepel", "pheatmap",
               "matrixStats", "data.table", "tidyr", "purrr", "igraph", "ggraph", "enrichR")
bioc_pkgs <- c("limma", "fgsea", "AnnotationDbi", "org.Dm.eg.db", "GO.db")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(limma)
  library(pheatmap)
  library(matrixStats)
  library(data.table)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  library(GO.db)
  library(enrichR)
  library(tidyr)
  library(igraph)
  library(ggraph)
  library(purrr)
})

select <- dplyr::select
filter <- dplyr::filter

# PATHS
BASE_DIR <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics"
TMT_DIR <- file.path(BASE_DIR, "TMT_all_from_psm_pdlike_ms3fix")
META_FILE <- file.path(BASE_DIR, "a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt")

# Output dirs
OUT_DIR <- file.path(BASE_DIR, "RESULTS_OSD514_ms3fix")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(OUT_DIR, "figs")
TAB_DIR <- file.path(OUT_DIR, "tables")

# LOAD METADATA
meta <- read_tsv(META_FILE, show_col_types = FALSE) %>%
  transmute(
    sample = str_replace_all(`Sample Name`, " ", "_"),
    tmt_run = factor(`Parameter Value[Run Number]`),
    condition = case_when(
      str_detect(`Sample Name`, "^Earth") ~ "Earth",
      str_detect(`Sample Name`, "^SF1g") ~ "SF1g",
      str_detect(`Sample Name`, "^SFug") ~ "SFug",
      TRUE ~ NA_character_
    ),
    sex = case_when(
      str_detect(`Sample Name`, "_M") ~ "Male",
      str_detect(`Sample Name`, "_F") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition), !is.na(sex)) %>%
  distinct(sample, .keep_all = TRUE)

meta$condition <- factor(meta$condition, levels = c("Earth", "SF1g", "SFug"))

# BUILD EXPRESSION MATRIX FROM MS3FIX TXT FILES
files <- list.files(TMT_DIR, pattern = "\\.txt$", full.names = TRUE)

read_tmt <- function(f) {
  sample_name <- basename(f) %>%
    str_remove("\\.txt$") %>%
    str_remove("_TMT[a-c]$")

  run_id <- basename(f) %>% str_extract("TMT[a-c](?=\\.txt$)")

  if (str_detect(sample_name, "^pool")) return(NULL)

  df <- read.delim(f, check.names = FALSE)

  ab_col <- grep("^Abundance:", colnames(df), value = TRUE)
  if (length(ab_col) != 1) {
    stop("Unexpected abundance columns in: ", basename(f))
  }

  df <- df %>%
    mutate(`# Unique Peptides` = suppressWarnings(as.numeric(`# Unique Peptides`))) %>%
    filter(`# Unique Peptides` >= 2)

  tibble(
    protein = df$Accession,
    sample = sample_name,
    run = run_id,
    value = suppressWarnings(as.numeric(df[[ab_col]]))
  ) %>%
    filter(!is.na(value), is.finite(value), value > 0) %>%
    group_by(protein, sample, run) %>%
    summarise(value = median(value, na.rm = TRUE), .groups = "drop")
}

expr_long <- map_dfr(files, read_tmt)

cat("Metadata samples:", nrow(meta), "\n")
cat("Total txt files:", length(files), "\n")
cat("Rows after read/filter:", nrow(expr_long), "\n")

# Keep proteins found across all TMT runs for robust run normalization
proteins_all_runs <- expr_long %>%
  group_by(protein) %>%
  summarise(n_runs = n_distinct(run), .groups = "drop") %>%
  filter(n_runs == 3) %>%
  pull(protein)

expr_long <- expr_long %>%
  filter(protein %in% proteins_all_runs)

# Within-run scaling
sample_medians <- expr_long %>%
  group_by(run, sample) %>%
  summarise(sample_median = median(value, na.rm = TRUE), .groups = "drop")

run_reference <- sample_medians %>%
  group_by(run) %>%
  summarise(run_median = median(sample_median, na.rm = TRUE), .groups = "drop")

scaling_tbl <- sample_medians %>%
  inner_join(run_reference, by = "run") %>%
  mutate(scale_factor = ifelse(sample_median > 0, run_median / sample_median, 1)) %>%
  select(run, sample, scale_factor)

expr_long <- expr_long %>%
  inner_join(scaling_tbl, by = c("run", "sample")) %>%
  mutate(value = value * scale_factor) %>%
  select(protein, sample, value)

expr_mat <- expr_long %>%
  pivot_wider(names_from = sample, values_from = value) %>%
  column_to_rownames("protein") %>%
  as.matrix()

common_samples <- intersect(colnames(expr_mat), meta$sample)

expr_limma <- expr_mat[, common_samples, drop = FALSE]
meta_limma <- meta %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

stopifnot(all(colnames(expr_limma) == meta_limma$sample))

expr_limma <- expr_limma[complete.cases(expr_limma), , drop = FALSE]
expr_limma <- log2(expr_limma + 1)
mode(expr_limma) <- "numeric"

write.csv(
  data.frame(protein_id = rownames(expr_limma), expr_limma, check.names = FALSE),
  file.path(TAB_DIR, "TMT_expression_matrix_ms3fix.csv"),
  row.names = FALSE
)

write.csv(
  meta_limma,
  file.path(TAB_DIR, "sample_metadata_ms3fix.csv"),
  row.names = FALSE
)

cat("Expression matrix:", nrow(expr_limma), "proteins x", ncol(expr_limma), "samples\n")

meta_limma$tmt_run <- factor(meta_limma$tmt_run)
meta_limma$sex <- factor(meta_limma$sex, levels = c("Female", "Male"))
meta_limma$group <- factor(
  paste(meta_limma$condition, meta_limma$sex, sep = "_"),
  levels = c("Earth_Female", "Earth_Male", "SF1g_Female", "SF1g_Male", "SFug_Female", "SFug_Male")
)

# PCA (raw and batch-corrected)
plot_pca <- function(mat, meta_df, title, subtitle, out_file) {
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  percent_var <- (pca$sdev^2) / sum(pca$sdev^2)

  pca_df <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(meta_df %>% select(sample, condition, sex, tmt_run), by = "sample")

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = tmt_run)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_bw(base_size = 12) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0("PC1 (", round(percent_var[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(percent_var[2] * 100, 1), "%)")
    )

  ggsave(out_file, p, width = 8, height = 6, dpi = 300)
}

plot_pca(
  mat = expr_limma,
  meta_df = meta_limma,
  title = "TMT proteomics PCA (ms3fix, raw)",
  subtitle = "Color = condition, shape = TMT run batch",
  out_file = file.path(FIG_DIR, "PCA_TMT_condition_batch_ms3fix_raw.png")
)

design_keep <- model.matrix(~ condition + sex, data = meta_limma)
expr_batch_corrected <- removeBatchEffect(
  expr_limma,
  batch = meta_limma$tmt_run,
  design = design_keep
)

write.csv(
  data.frame(protein_id = rownames(expr_batch_corrected), expr_batch_corrected, check.names = FALSE),
  file.path(TAB_DIR, "TMT_expression_matrix_ms3fix_batch_corrected.csv"),
  row.names = FALSE
)

plot_pca(
  mat = expr_batch_corrected,
  meta_df = meta_limma,
  title = "TMT proteomics PCA (ms3fix, batch-corrected)",
  subtitle = "Batch removed with limma::removeBatchEffect (batch = tmt_run)",
  out_file = file.path(FIG_DIR, "PCA_TMT_condition_batch_ms3fix_batch_corrected.png")
)

# LIMMA model specs
CONTRAST_NAMES <- c("SF1g_vs_Earth", "SFug_vs_Earth", "SF1g_vs_SFug")

LIMMA_PADJ_CUT <- 0.05
LIMMA_LOGFC_CUT <- 0.5
LIMMA_FC_CUT <- 2^LIMMA_LOGFC_CUT

MODEL_SPECS <- tibble::tribble(
  ~model_id,              ~formula_str,
  "condition_only",      "~0 + condition",
  "condition_batch",     "~0 + condition + tmt_run",
  "condition_batch_sex", "~0 + condition + tmt_run + sex"
)

MODEL_OUT_DIR <- file.path(TAB_DIR, "model_outputs")
MODEL_FIG_DIR <- file.path(FIG_DIR, "model_outputs")
dir.create(MODEL_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

build_design <- function(formula_str, meta_df, condition_levels = levels(meta_df$condition)) {
  design <- model.matrix(as.formula(formula_str), data = meta_df)
  cond_cols <- grep("^condition", colnames(design))

  if (length(cond_cols) != length(condition_levels)) {
    stop(
      "Could not find expected number of condition columns (",
      length(condition_levels), ") in design for formula: ",
      formula_str
    )
  }

  colnames(design)[cond_cols] <- condition_levels
  design
}

pairwise_conditions_from_contrast <- function(contrast_name) {
  parts <- strsplit(contrast_name, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop("Unexpected contrast format: ", contrast_name)
  parts
}

save_volcano <- function(fit_obj, contrast_name,
                         fig_dir,
                         title_prefix = "",
                         pval_cut = LIMMA_PADJ_CUT,
                         logfc_cut = LIMMA_LOGFC_CUT,
                         top_n = 5) {

  tt <- topTable(fit_obj, coef = contrast_name,
                 number = Inf, adjust.method = "BH") %>%
    rownames_to_column("protein_id") %>%
    filter(!is.na(P.Value), !is.na(adj.P.Val), !is.na(logFC)) %>%
    mutate(
      sig = case_when(
        adj.P.Val < pval_cut & logFC >= logfc_cut  ~ "Up",
        adj.P.Val < pval_cut & logFC <= -logfc_cut ~ "Down",
        TRUE ~ "NotSig"
      ),
      sig = factor(sig, levels = c("Down", "NotSig", "Up"))
    )

  top_labels <- bind_rows(
    tt %>% filter(sig == "Up") %>% arrange(desc(logFC)) %>% slice_head(n = top_n),
    tt %>% filter(sig == "Down") %>% arrange(logFC) %>% slice_head(n = top_n)
  )

  p <- ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
    geom_text_repel(
      data = top_labels,
      aes(label = protein_id),
      max.overlaps = 50,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(Down = "blue", NotSig = "grey70", Up = "red"), drop = FALSE) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0("Volcano Plot: ", title_prefix, contrast_name),
      x = "Log2 Fold Change",
      y = expression(-log[10](adjusted~italic(p)))
    )

  ggsave(file.path(fig_dir, paste0("Volcano_", contrast_name, ".png")),
         p, width = 8, height = 6, dpi = 300)
}

run_limma_model <- function(model_id, formula_str) {
  cat("\nRunning limma model:", model_id, "with", formula_str, "\n")

  model_tab_dir <- file.path(MODEL_OUT_DIR, model_id)
  model_fig_dir <- file.path(MODEL_FIG_DIR, model_id)
  dir.create(model_tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(model_fig_dir, recursive = TRUE, showWarnings = FALSE)

  results <- list()
  sig <- list()

  for (cn in CONTRAST_NAMES) {
    cond_pair <- pairwise_conditions_from_contrast(cn)

    meta_sub <- meta_limma %>%
      filter(condition %in% cond_pair) %>%
      droplevels()

    meta_sub$condition <- factor(meta_sub$condition, levels = cond_pair)

    expr_sub <- expr_limma[, meta_sub$sample, drop = FALSE]

    design <- build_design(
      formula_str = formula_str,
      meta_df = meta_sub,
      condition_levels = cond_pair
    )

    contrast_matrix <- makeContrasts(
      contrasts = paste0(cond_pair[[1]], " - ", cond_pair[[2]]),
      levels = design
    )
    colnames(contrast_matrix) <- cn

    fit <- lmFit(expr_sub, design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)

    write.csv(
      data.frame(sample = meta_sub$sample, design, check.names = FALSE),
      file.path(model_tab_dir, paste0("design_matrix_", cn, ".csv")),
      row.names = FALSE
    )

    write.csv(
      meta_sub %>% dplyr::select(sample, condition, sex, tmt_run),
      file.path(model_tab_dir, paste0("model_metadata_", cn, ".csv")),
      row.names = FALSE
    )

    tt <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH") %>%
      rownames_to_column("protein_id")

    tt_sig <- tt %>% filter(adj.P.Val < LIMMA_PADJ_CUT, abs(logFC) >= LIMMA_LOGFC_CUT)

    results[[cn]] <- tt
    sig[[cn]] <- tt_sig

    write.csv(
      tt,
      file.path(model_tab_dir, paste0("Limma_", cn, "_results.csv")),
      row.names = FALSE
    )

    write.csv(
      tt_sig,
      file.path(model_tab_dir, paste0("Limma_", cn, "_significant.csv")),
      row.names = FALSE
    )

    save_volcano(
      fit_obj = fit2,
      contrast_name = cn,
      fig_dir = model_fig_dir,
      title_prefix = paste0(model_id, " | ")
    )
  }

  summary_tbl <- tibble(
    model = model_id,
    contrast = CONTRAST_NAMES,
    n_significant = vapply(sig[CONTRAST_NAMES], nrow, integer(1))
  )

  write.csv(
    summary_tbl,
    file.path(model_tab_dir, "significant_counts_summary.csv"),
    row.names = FALSE
  )

  summary_tbl
}

run_sex_stratified_limma <- function() {
  model_id <- "sex_stratified_condition_batch"
  cat("\nRunning limma model:", model_id, "with ~0 + group + tmt_run\n")

  model_tab_dir <- file.path(MODEL_OUT_DIR, model_id)
  model_fig_dir <- file.path(MODEL_FIG_DIR, model_id)
  dir.create(model_tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(model_fig_dir, recursive = TRUE, showWarnings = FALSE)

  sex_contrast_groups <- list(
    SFug_females_vs_Earth_females = c("SFug_Female", "Earth_Female"),
    SFug_males_vs_Earth_males = c("SFug_Male", "Earth_Male"),
    SF1g_females_vs_Earth_females = c("SF1g_Female", "Earth_Female"),
    SF1g_males_vs_Earth_males = c("SF1g_Male", "Earth_Male")
  )

  sex_contrasts <- names(sex_contrast_groups)
  summary_tbl <- tibble(
    model = character(),
    contrast = character(),
    n_significant_fdr05 = integer(),
    n_significant_fdr05_lfc1 = integer()
  )

  summary_detailed_tbl <- tibble(
    model = character(),
    contrast = character(),
    n_tested = integer(),
    n_significant_fdr05 = integer(),
    n_significant_fdr05_lfc1 = integer(),
    n_up_fdr05_lfc1 = integer(),
    n_down_fdr05_lfc1 = integer()
  )

  for (cn in sex_contrasts) {
    group_pair <- sex_contrast_groups[[cn]]

    meta_sub <- meta_limma %>%
      filter(group %in% group_pair) %>%
      droplevels()

    meta_sub$group <- factor(meta_sub$group, levels = group_pair)

    expr_sub <- expr_limma[, meta_sub$sample, drop = FALSE]

    design <- model.matrix(~0 + group + tmt_run, data = meta_sub)
    group_cols <- grep("^group", colnames(design))
    colnames(design)[group_cols] <- group_pair

    contrast_matrix <- makeContrasts(
      contrasts = paste0(group_pair[[1]], " - ", group_pair[[2]]),
      levels = design
    )
    colnames(contrast_matrix) <- cn

    fit <- lmFit(expr_sub, design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)

    write.csv(
      data.frame(sample = meta_sub$sample, design, check.names = FALSE),
      file.path(model_tab_dir, paste0("design_matrix_", cn, ".csv")),
      row.names = FALSE
    )

    write.csv(
      meta_sub %>% dplyr::select(sample, condition, sex, tmt_run, group),
      file.path(model_tab_dir, paste0("model_metadata_", cn, ".csv")),
      row.names = FALSE
    )

    tt <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH") %>%
      rownames_to_column("protein_id")

    tt_sig_fdr <- tt %>% filter(adj.P.Val < LIMMA_PADJ_CUT)
    tt_sig <- tt_sig_fdr %>% filter(abs(logFC) >= LIMMA_LOGFC_CUT)

    write.csv(
      tt,
      file.path(model_tab_dir, paste0("Limma_", cn, "_results.csv")),
      row.names = FALSE
    )

    write.csv(
      tt_sig_fdr,
      file.path(model_tab_dir, paste0("Limma_", cn, "_significant_fdr_only.csv")),
      row.names = FALSE
    )

    write.csv(
      tt_sig,
      file.path(model_tab_dir, paste0("Limma_", cn, "_significant.csv")),
      row.names = FALSE
    )

    save_volcano(
      fit_obj = fit2,
      contrast_name = cn,
      fig_dir = model_fig_dir,
      title_prefix = paste0(model_id, " | ")
    )

    summary_tbl <- bind_rows(
      summary_tbl,
      tibble(
        model = model_id,
        contrast = cn,
        n_significant_fdr05 = nrow(tt_sig_fdr),
        n_significant_fdr05_lfc1 = nrow(tt_sig)
      )
    )

    summary_detailed_tbl <- bind_rows(
      summary_detailed_tbl,
      tibble(
        model = model_id,
        contrast = cn,
        n_tested = nrow(tt),
        n_significant_fdr05 = nrow(tt_sig_fdr),
        n_significant_fdr05_lfc1 = nrow(tt_sig),
        n_up_fdr05_lfc1 = sum(tt_sig$logFC > 0, na.rm = TRUE),
        n_down_fdr05_lfc1 = sum(tt_sig$logFC < 0, na.rm = TRUE)
      )
    )
  }

  write.csv(
    summary_tbl,
    file.path(model_tab_dir, "significant_counts_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    summary_detailed_tbl,
    file.path(model_tab_dir, "significant_counts_direction_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    summary_tbl,
    file.path(TAB_DIR, "sex_stratified_significant_counts_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    summary_detailed_tbl,
    file.path(TAB_DIR, "sex_stratified_significant_counts_direction_summary.csv"),
    row.names = FALSE
  )

  summary_tbl
}

model_summaries <- list()
for (i in seq_len(nrow(MODEL_SPECS))) {
  model_id <- MODEL_SPECS$model_id[i]
  formula_str <- MODEL_SPECS$formula_str[i]
  model_summaries[[model_id]] <- run_limma_model(model_id, formula_str)
}

all_model_summary <- bind_rows(model_summaries)
write.csv(
  all_model_summary,
  file.path(MODEL_OUT_DIR, "all_model_significant_counts_summary.csv"),
  row.names = FALSE
)

sex_model_summary <- run_sex_stratified_limma()
write.csv(
  sex_model_summary,
  file.path(MODEL_OUT_DIR, "sex_stratified_significant_counts_summary.csv"),
  row.names = FALSE
)

PRIMARY_MODEL <- "condition_batch_sex"
primary_tab_dir <- file.path(MODEL_OUT_DIR, PRIMARY_MODEL)

results_SF1g_vs_Earth <- read.csv(
  file.path(primary_tab_dir, "Limma_SF1g_vs_Earth_results.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

results_SFug_vs_Earth <- read.csv(
  file.path(primary_tab_dir, "Limma_SFug_vs_Earth_results.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

results_SF1g_vs_SFug <- read.csv(
  file.path(primary_tab_dir, "Limma_SF1g_vs_SFug_results.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

top_proteins <- results_SF1g_vs_Earth %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(protein_id)

top_proteins <- top_proteins[top_proteins %in% rownames(expr_limma)]
heat_mat <- expr_limma[top_proteins, , drop = FALSE]

ann <- meta_limma %>%
  dplyr::select(sample, condition, sex, tmt_run) %>%
  column_to_rownames("sample")

ord <- order(ann$condition, ann$sex, ann$tmt_run)
heat_mat <- heat_mat[, ord, drop = FALSE]
ann <- ann[ord, , drop = FALSE]

heat_z <- t(scale(t(heat_mat)))
heat_z[is.na(heat_z)] <- 0

png(file.path(FIG_DIR, "Heatmap_Top50_ms3fix.png"), 2200, 1500, res = 200)
pheatmap(
  heat_z,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann,
  show_rownames = TRUE,
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA,
  main = "Top 50 Differential Proteins across Earth, SF1g, SFug (ms3fix)"
)
dev.off()

write.csv(results_SF1g_vs_Earth, file.path(TAB_DIR, "Limma_SF1g_vs_Earth_results.csv"), row.names = FALSE)
write.csv(results_SFug_vs_Earth, file.path(TAB_DIR, "Limma_SFug_vs_Earth_results.csv"), row.names = FALSE)
write.csv(results_SF1g_vs_SFug, file.path(TAB_DIR, "Limma_SF1g_vs_SFug_results.csv"), row.names = FALSE)

sig_SF1g_vs_Earth <- results_SF1g_vs_Earth %>% filter(adj.P.Val < 0.05, abs(logFC) >= 1)
sig_SFug_vs_Earth <- results_SFug_vs_Earth %>% filter(adj.P.Val < 0.05, abs(logFC) >= 1)
sig_SF1g_vs_SFug <- results_SF1g_vs_SFug %>% filter(adj.P.Val < 0.05, abs(logFC) >= 1)

write.csv(sig_SF1g_vs_Earth, file.path(TAB_DIR, "Limma_SF1g_vs_Earth_significant.csv"), row.names = FALSE)
write.csv(sig_SFug_vs_Earth, file.path(TAB_DIR, "Limma_SFug_vs_Earth_significant.csv"), row.names = FALSE)
write.csv(sig_SF1g_vs_SFug, file.path(TAB_DIR, "Limma_SF1g_vs_SFug_significant.csv"), row.names = FALSE)

cat("SF1g vs Earth:", nrow(sig_SF1g_vs_Earth), "significant proteins\n")
cat("SFug vs Earth:", nrow(sig_SFug_vs_Earth), "significant proteins\n")
cat("SF1g vs SFug:", nrow(sig_SF1g_vs_SFug), "significant proteins\n")

## GSEA

GSEA_DIR <- file.path(OUT_DIR, "GSEA")
TBL_DIR <- TAB_DIR

dir.create(file.path(GSEA_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(GSEA_DIR, "figs"), recursive = TRUE, showWarnings = FALSE)

PADJ_CUT <- 0.05
LFC_CUT <- 1
CONTRASTS <- c("SF1g_vs_Earth", "SFug_vs_Earth", "SF1g_vs_SFug")

`%||%` <- function(a, b) if (!is.null(a) && length(a)) a else b
wrap_terms <- function(x, w = 40) stringr::str_wrap(x, w)

load_de <- function(tag) {
  f <- file.path(TBL_DIR, paste0("Limma_", tag, "_results.csv"))
  if (!file.exists(f)) stop("Missing DE file: ", f)

  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(df) <- df$protein_id
  df
}

build_ranks <- function(res) {
  ids <- rownames(res)
  pv <- res$P.Value %||% res$adj.P.Val
  lfc <- res$logFC %||% 0

  pv[!is.finite(pv)] <- 1
  lfc[!is.finite(lfc)] <- 0

  score <- sign(lfc) * (-log10(pmax(pv, 1e-300)))

  keep <- is.finite(score) & nzchar(ids)
  ranks <- tapply(score[keep], ids[keep], mean)

  sort(ranks, decreasing = TRUE)
}

build_go_bp_sets <- function(ranks) {
  genes <- names(ranks)

  m <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = genes,
    keytype = "UNIPROT",
    columns = c("GOALL", "ONTOLOGYALL")
  )

  m <- m[m$ONTOLOGYALL == "BP", ]
  if (nrow(m) == 0) return(NULL)

  paths <- split(m$UNIPROT, m$GOALL)
  paths <- lapply(paths, function(x) unique(na.omit(x)))

  lens <- lengths(paths)
  paths[lens >= 10 & lens <= 500]
}

get_go_terms <- function(ids) {
  if (!length(ids)) return(setNames(character(0), character(0)))

  tbl <- AnnotationDbi::select(
    GO.db,
    keys = ids,
    keytype = "GOID",
    columns = "TERM"
  )

  tbl <- unique(tbl)
  setNames(tbl$TERM, tbl$GOID)
}

get_sig <- function(res) {
  res <- res[
    !is.na(res$adj.P.Val) &
      res$adj.P.Val < PADJ_CUT &
      !is.na(res$logFC) &
      abs(res$logFC) >= LFC_CUT,
  ]
  unique(rownames(res))
}

uniprot_to_symbol <- function(ids) {
  tbl <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = ids,
    keytype = "UNIPROT",
    columns = "SYMBOL"
  )

  tbl <- tbl[!duplicated(tbl$UNIPROT), ]
  map <- setNames(tbl$SYMBOL, tbl$UNIPROT)

  out <- map[ids]
  out[is.na(out)] <- ids[is.na(out)]
  out
}

save_fgsea <- function(dt, tag) {
  if (is.null(dt) || nrow(dt) == 0) return()

  dt <- dt[order(padj)][1:min(15, nrow(dt))]
  term_map <- get_go_terms(dt$pathway)

  dt$term <- term_map[dt$pathway]
  dt$term[is.na(dt$term)] <- dt$pathway
  dt$label <- wrap_terms(dt$term)

  p <- ggplot(dt, aes(reorder(label, NES), NES, fill = -log10(padj))) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("GSEA:", tag))

  ggsave(file.path(GSEA_DIR, "figs", paste0("fgsea_", tag, ".png")),
         p, width = 10, height = 6)
}

save_enrichment_plots <- function(pathways, ranks, fg_dt, tag) {
  if (is.null(fg_dt) || nrow(fg_dt) == 0) return()

  top <- fg_dt[order(padj)][1:min(5, nrow(fg_dt))]
  term_map <- get_go_terms(top$pathway)

  for (i in seq_len(nrow(top))) {
    pw <- top$pathway[i]
    term_name <- term_map[pw]
    if (is.na(term_name)) term_name <- pw

    p <- plotEnrichment(pathways[[pw]], ranks) +
      labs(title = paste0(tag, " | ", term_name))

    ggsave(
      file.path(GSEA_DIR, "figs", paste0("enrichment_", tag, "_", i, ".png")),
      p,
      width = 8,
      height = 6
    )
  }
}

run_enrichr <- function(tag, genes) {
  if (!length(genes)) return(NULL)

  setEnrichrSite("FlyEnrichr")

  dbs <- grep("^GO_", listEnrichrDbs()$libraryName, value = TRUE)
  enrich <- enrichr(genes, dbs)

  for (d in names(enrich)) {
    fwrite(as.data.table(enrich[[d]]),
           file.path(GSEA_DIR, "tables", paste0(tag, "_", d, ".csv")))
  }

  enrich
}

run_network <- function(enrich_bp, tag) {
  if (is.null(enrich_bp) || nrow(enrich_bp) == 0) return()

  top <- enrich_bp %>%
    arrange(Adjusted.P.value) %>%
    head(12)

  edge_list <- lapply(seq_len(nrow(top)), function(i) {
    genes <- unlist(strsplit(top$Genes[i], ";"))
    genes <- genes[genes != "" & !is.na(genes)]

    if (!length(genes)) return(NULL)

    data.frame(
      from = rep(top$Term[i], length(genes)),
      to = genes,
      stringsAsFactors = FALSE
    )
  })

  edge_list <- edge_list[!sapply(edge_list, is.null)]
  if (!length(edge_list)) return()

  edges <- do.call(rbind, edge_list)
  g <- graph_from_data_frame(edges)

  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.3) +
    geom_node_point(size = 3) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    theme_void()

  ggsave(file.path(GSEA_DIR, "figs", paste0("network_", tag, ".png")),
         p, width = 10, height = 8)
}

run_analysis <- function(tag) {
  message("\nRUNNING: ", tag)

  res <- load_de(tag)
  ranks_uniprot <- build_ranks(res)

  pathways <- build_go_bp_sets(ranks_uniprot)
  if (is.null(pathways) || length(pathways) < 10) {
    message("Skipping ", tag, ": too few pathways")
    return(NULL)
  }

  fg <- fgseaMultilevel(pathways, ranks_uniprot)
  fg_dt <- as.data.table(fg)

  if (nrow(fg_dt) > 0) {
    fwrite(fg_dt, file.path(GSEA_DIR, "tables", paste0("fgsea_", tag, ".csv")))
    save_fgsea(fg_dt, tag)
    save_enrichment_plots(pathways, ranks_uniprot, fg_dt, tag)
  }

  sig <- get_sig(res)
  sig_symbols <- uniprot_to_symbol(sig)

  enrich <- run_enrichr(tag, sig_symbols)

  if (!is.null(enrich[["GO_Biological_Process_2018"]])) {
    run_network(enrich[["GO_Biological_Process_2018"]], tag)
  }
}

for (t in CONTRASTS) {
  tryCatch(run_analysis(t), error = function(e) message("ERROR: ", e$message))
}

cat("\nDONE → ", OUT_DIR, "\n")
