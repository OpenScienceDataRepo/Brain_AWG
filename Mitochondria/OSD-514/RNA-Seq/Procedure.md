# Extract Raw Counts

OSD-514 did not have a clean raw counts file to run in DESeq2. It had an "unnormalized counts" file but that had non-integer values, meaning it was a convenience matrix of RSEM “expected counts” (decimals).

The 24 files compiled are also not raw counts files. The *.genes.results files from RSEM contain expected counts (fractional). We used them because DESeq2 supports this workflow via tximport. We gave DESeq2 the expected counts plus the average transcript lengths that tximport extracts. DESeq2 then uses these as offsets (length-aware normalization).

Just rounding the “Unnormalized Counts” CSV throws away the information RSEM provides. It biases low counts, meaning small fractional differences can flip to 0 or 1 after rounding, distorting dispersion and p-values. It also breaks length-aware normalization since tximport+DESeq2 can properly correct for gene length and library size.

Because of this, the 24 raw count results were downloaded for each sample in each condition and collapsed into one matrix manually using the below code.

```
# Set Working Directory first

# Set the output folder
out_dir <- "Collapsed_Counts"
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# Find the 24 gene files
files <- list.files("RawCounts", pattern="\\.genes\\.results$", full.names=TRUE)
stopifnot(length(files) == 24)

# Helper to read one file -> named numeric vector of expected_count
read_one <- function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  # RSEM gene files have 'gene_id' and 'expected_count' columns
  if (!all(c("gene_id","expected_count") %in% colnames(d))) {
    stop("File lacks required columns: ", basename(f))
  }
  v <- as.numeric(d$expected_count)
  names(v) <- d$gene_id
  v
}

# Read all samples
lst <- lapply(files, read_one)

# Make matrix (fill missing with 0)
all_genes <- Reduce(union, lapply(lst, names))
mat <- sapply(lst, function(v) { x <- v[all_genes]; x[is.na(x)] <- 0; x })
rownames(mat) <- all_genes

# Nice column names from filenames
samp <- sub("\\.genes\\.results$", "", basename(files))
colnames(mat) <- samp

# Write the decimal "expected counts" matrix (for record/tximport)
write.csv(mat,
          file.path(out_dir, "tables", "OSD514_RSEM_expected_counts.csv"),
          quote = FALSE)

# Also write a rounded integer matrix that DESeq2 will accept directly
mat_int <- round(mat)
storage.mode(mat_int) <- "integer"
write.csv(mat_int,
          file.path(out_dir, "tables", "counts_DESeq_ready.csv"),
          quote = FALSE)

# Make a tiny metadata from filenames (condition_group, sex)
grp_tok <- sub("^GLDS-514_rna-seq_([^_]+)_.*$", "\\1", samp)       # Earth / SF1g / SFug
sex_tok <- sub("^GLDS-514_rna-seq_[^_]+_([MF])\\d+_.*$", "\\1", samp) # M/F
cond_map <- c(Earth="EARTH", SF1g="SPACEFLIGHT_1G", SFug="SPACEFLIGHT_MICROGRAVITY")
meta <- data.frame(
  sample = samp,
  condition_group = factor(cond_map[grp_tok],
                           levels = c("EARTH","SPACEFLIGHT_1G","SPACEFLIGHT_MICROGRAVITY")),
  sex = factor(ifelse(sex_tok=="F","FEMALE","MALE"))
)
rownames(meta) <- meta$sample
write.csv(meta, file.path(out_dir, "tables", "metadata_from_filenames.csv"),
          row.names = TRUE, quote = FALSE)

cat("Wrote:\n",
    "- tables/OSD514_RSEM_expected_counts.csv\n",
    "- tables/counts_DESeq_ready.csv (rounded; plug into DESeq2)\n",
    "- tables/metadata_from_filenames.csv (optional)\n", sep = "")
```

# Quality Control
```
# 0) Working directory
setwd("~/Desktop/ADBR Mito/OSD-514")

# 1) Load packages
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tibble)
  library(stringr)
})

# 2) Output dirs
OUT_DIR <- "RESULTS_OSD514"
dir.create(file.path(OUT_DIR, "figs"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

# 3) Find the gene result files
files <- list.files(path=".", pattern="\\.genes\\.results$", full.names=TRUE, recursive=TRUE)
if (length(files) == 0) {
  stop("No *.genes.results files found. Make sure you've downloaded the 24 '...genes.results' files from OSDR.")
}
message("Found ", length(files), " *.genes.results files.")
# If you want to enforce all 24 samples, uncomment:
# if (length(files) != 24) stop("Expected 24 files; found ", length(files), ".")

# 4) Parse tokens from filenames to build metadata
b <- basename(files)
# GLDS-514_rna-seq_<SFug|SF1g|Earth>_<M|F><rep>_..._(CRRA#######)-..._(HV... )_L<lane>.genes.results
m <- stringr::str_match(b, "rna-seq_(SFug|SF1g|Earth)_([MF])(\\d).*_(CRRA\\d+)\\-[^_]+_(HV\\w+)_L(\\d)")
if (any(is.na(m))) {
  bad <- b[apply(is.na(m), 1, any)]
  stop("Filename parsing failed for these files:\n  - ", paste(bad, collapse="\n  - "),
       "\nIf names differ, adjust the regex in str_match().")
}

cond_map <- c(SFug="SPACEFLIGHT_MICROGRAVITY", SF1g="SPACEFLIGHT_1G", Earth="EARTH")
sex_map  <- c(M="MALE", F="FEMALE")

condition_group <- cond_map[m[,2]]
sex             <- sex_map[m[,3]]
replicate       <- as.integer(m[,4])
crra            <- m[,5]
flowcell        <- m[,6]
lane            <- paste0("L", m[,7])

# Human-readable sample names
samp <- paste(condition_group, sex, replicate, crra, flowcell, lane)
samp <- make.unique(samp)
names(files) <- samp

meta <- data.frame(
  condition_group = factor(condition_group, levels=c("EARTH","SPACEFLIGHT_1G","SPACEFLIGHT_MICROGRAVITY")),
  sex             = factor(sex, levels=c("FEMALE","MALE")),
  replicate       = replicate,
  crra            = crra,
  flowcell        = flowcell,
  lane            = lane,
  row.names       = samp,
  check.names     = FALSE
)
write.csv(meta, file.path(OUT_DIR,"tables","sample_table.csv"))

# 5) Import RSEM with tximport
# Unscaled expected counts for a "library size" proxy (QC only)
txi_raw <- tximport(files, type = "rsem", countsFromAbundance = "no")
lib <- colSums(txi_raw$counts)
write.csv(data.frame(sample = names(lib), libsize_sum_expected_counts = lib),
          file.path(OUT_DIR,"tables","01_library_sizes.csv"), row.names = FALSE)

png(file.path(OUT_DIR,"figs","01_library_sizes_boxplot.png"), width=1200, height=900, res=150)
boxplot(log10(lib), ylab = "log10 sum(expected counts)", main = "Library size proxy (log10)")
dev.off()

flag_low <- names(lib)[lib < median(lib)/3]
if (length(flag_low)) message("Flagged low-depth (~≥3x below median): ", paste(flag_low, collapse = ", "))

# Gene-level expected counts for DESeq2 (no rounding)
txi <- tximport(files, type = "rsem", countsFromAbundance = "no")

# 6) Low-count filter (≥10 in ≥20% samples)
keep <- rowSums(txi$counts >= 10) >= ceiling(ncol(txi$counts) * 0.2)
txi$counts    <- txi$counts[keep, , drop = FALSE]
txi$abundance <- txi$abundance[keep, , drop = FALSE]
txi$length    <- txi$length[keep, , drop = FALSE]

## ---- PRINT "02_filtering_summary.csv" IN CONSOLE ----
filtering_summary <- data.frame(
  kept_genes = nrow(txi$counts),
  rule       = ">=10 in >=20% samples",
  check.names = FALSE
)
cat("\n===== 02_filtering_summary.csv (printed) =====\n")
print(filtering_summary, row.names = FALSE)

# 7) Build DESeq2 object + normalization
design_formula <- if (any(is.na(meta$sex))) ~ condition_group else ~ sex + condition_group
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = design_formula)
dds <- estimateSizeFactors(dds)

sf <- tryCatch(sizeFactors(dds), error = function(e) NULL)
if (!is.null(sf) && length(sf) == ncol(dds)) {
  write.csv(data.frame(sample = colnames(dds), sizeFactor = sf),
            file.path(OUT_DIR,"tables","03_size_factors.csv"), row.names = FALSE)
} else {
  writeLines(
    "Note: DESeq2 is using normalization factors (offsets) rather than simple sizeFactors; sizeFactors() is NULL/unused in this setting.",
    con = file.path(OUT_DIR,"tables","03_normalization_note.txt")
  )
}

# 8) VST → PCA  (with 95% ellipses per condition)
vsd <- vst(dds, blind = TRUE)
pca_df <- plotPCA(vsd, intgroup = c("condition_group","sex"), returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2)) +
  # translucent ellipses for each condition
  stat_ellipse(aes(fill = condition_group, group = condition_group),
               type = "norm", level = 0.95, geom = "polygon",
               alpha = 0.15, color = NA) +
  # points on top
  geom_point(aes(color = condition_group, shape = sex), size = 3) +
  labs(x = paste0("PC1 (", percentVar[1], "%)"),
       y = paste0("PC2 (", percentVar[2], "%)"),
       title = paste0("PCA on VST data (design: ", deparse(design_formula), ")")) +
  theme_classic()

ggsave(file.path(OUT_DIR,"figs","04_pca.png"), p_pca, width=7, height=6, dpi=300, bg="white")

# 9) Sample-to-sample distance heatmap
#    Order *explicitly* by condition first, then sex, then replicate (NO gaps)
mat <- assay(vsd)

# Short labels for columns
map_cond <- c(EARTH="E", SPACEFLIGHT_1G="1G", SPACEFLIGHT_MICROGRAVITY="MG")
cond_code <- map_cond[as.character(meta$condition_group)]
sex_code  <- substr(as.character(meta$sex), 1, 1)
short_lab <- paste(cond_code, sex_code, meta$replicate, meta$lane, sep="_")

# Explicit ordinal index: condition → sex → replicate
cond_levels <- levels(meta$condition_group)    # E, 1G, MG
sex_levels  <- levels(meta$sex)                # FEMALE, MALE
ord <- unlist(lapply(cond_levels, function(cc) {
  idx_c <- which(meta$condition_group == cc)
  idx_c[order(factor(meta$sex[idx_c], levels = sex_levels),
              meta$replicate[idx_c])]
}), use.names = FALSE)

# Reorder matrices and annotations
D    <- as.matrix(dist(t(mat)))[ord, ord]
ann  <- meta[ord, c("condition_group","sex"), drop = FALSE]
lbls <- short_lab[ord]

# Save a key for labels
write.csv(data.frame(short=short_lab, sample=rownames(meta),
                     condition=meta$condition_group, sex=meta$sex,
                     crra=meta$crra, flowcell=meta$flowcell, lane=meta$lane),
          file.path(OUT_DIR,"tables","sample_label_key.csv"), row.names=FALSE)

png(file.path(OUT_DIR,"figs","05_sample_distance_heatmap.png"),
    width = 2000, height = 1600, res = 200)
pheatmap(D,
  annotation_col = ann,
  cluster_rows   = FALSE,   # keep our explicit order
  cluster_cols   = FALSE,   # keep our explicit order
  show_rownames  = FALSE,
  labels_col     = lbls,
  fontsize_col   = 10,
  border_color   = NA
)
dev.off()

# 10) Optional: batch token & summary
meta$batch <- factor(paste(meta$flowcell, meta$lane, sep="_"))
tab <- table(meta$condition_group, meta$batch)
write.csv(as.data.frame(tab), file.path(OUT_DIR,"tables","06_condition_by_batch.csv"), row.names=FALSE)

final_design <- if (any(colSums(tab > 0) >= 2) && nlevels(meta$batch) > 1) {
  "~ batch + sex + condition_group"
} else {
  deparse(design_formula)
}

## ---- PRINT "00_step2_summary.txt" IN CONSOLE ----
summary_lines <- c(
  sprintf("Library size proxy (sum expected counts): %.2fM–%.2fM", min(lib)/1e6, max(lib)/1e6),
  sprintf("Kept genes: %d", nrow(txi$counts)),
  sprintf("PCA: PC1=%d%%, PC2=%d%%", percentVar[1], percentVar[2]),
  paste("Design for DE:", final_design)
)
cat("\n===== 00_step2_summary.txt (printed) =====\n",
    paste(summary_lines, collapse = "\n"),
    "\n", sep = "")


```

# DESeq2
```
# Safety checks
stopifnot(exists("txi"), exists("meta"), exists("OUT_DIR"))
stopifnot(all(c("counts","abundance","length") %in% names(txi)))
stopifnot(all(colnames(txi$counts) == rownames(meta)))

# Ensure output folders exist
dir.create(file.path(OUT_DIR, "figs"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

# install shrinkers if missing
if (!requireNamespace("apeglm", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("apeglm", ask = FALSE, update = FALSE)
}

# 3.1 Design (reference + available covariates)
# Use EARTH as the control/reference
meta$condition_group <- stats::relevel(meta$condition_group, ref = "EARTH")
if (!is.factor(meta$sex))   meta$sex   <- factor(meta$sex)

# If batch not present, create from flowcell+lane (from Step 2)
if (!("batch" %in% colnames(meta)) && all(c("flowcell","lane") %in% colnames(meta))) {
  meta$batch <- factor(paste(meta$flowcell, meta$lane, sep = "_"))
}

use_batch <- "batch" %in% colnames(meta) && nlevels(droplevels(meta$batch)) > 1
use_sex   <- "sex"   %in% colnames(meta) && nlevels(droplevels(meta$sex))   > 1

if (use_batch && use_sex) {
  design_formula <- ~ batch + sex + condition_group
} else if (use_batch && !use_sex) {
  design_formula <- ~ batch + condition_group
} else if (!use_batch && use_sex) {
  design_formula <- ~ sex + condition_group
} else {
  design_formula <- ~ condition_group
}
message("Design: ", deparse(design_formula))

# Build & fit
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = design_formula)
dds <- DESeq(dds)  # size factors, dispersions, GLM

# Dispersion trend plot
png(file.path(OUT_DIR, "figs", "de_01_dispersion_trend.png"), width=1200, height=900, res=150)
plotDispEsts(dds)
dev.off()

# Peek at available coefficients
print(resultsNames(dds))

# Robust contrast runner: always saves unshrunk; shrinks via apeglm (relevel if needed), else optional ashr
run_contrast <- function(var, levelA, levelB, tag) {
  # 1) Unshrunk results (for padj)
  res <- results(dds, contrast = c(var, levelA, levelB), alpha = 0.05)
  res_df <- as.data.frame(res); res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$padj), ]
  out_unshr <- file.path(OUT_DIR, "tables", paste0("DE_", tag, "_unshrunk.csv"))
  write.csv(res_df, out_unshr, row.names = FALSE)

  # 2) Shrink LFCs for plotting/ranking
  coef_name <- paste0(var, "_", levelA, "_vs_", levelB)
  res_shr_df <- NULL

  # helper for MA plots
  ma_plot <- function(obj, label) {
    png(file.path(OUT_DIR, "figs", paste0("de_03_MA_", tag, "_", label, ".png")), width=1200, height=900, res=150)
    plotMA(obj, ylim = c(-4,4), main = paste("MA (shrunk:", label, "):", tag))
    dev.off()
  }

  if (coef_name %in% resultsNames(dds)) {
    # Direct coef -> apeglm
    res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    res_shr_df <- as.data.frame(res_shr); res_shr_df$gene <- rownames(res_shr_df)
    write.csv(res_shr_df, file.path(OUT_DIR, "tables", paste0("DE_", tag, "_LFCshrunk_apeglm.csv")), row.names = FALSE)
    ma_plot(res_shr, "apeglm")
  } else {
    # Relevel a TEMP copy so levelB becomes reference (turns contrast into a coef)
    dds_tmp <- dds
    colData(dds_tmp)[[var]] <- stats::relevel(colData(dds_tmp)[[var]], ref = levelB)
    dds_tmp <- DESeq(dds_tmp, quiet = TRUE)
    coef_name2 <- paste0(var, "_", levelA, "_vs_", levelB)
    if (coef_name2 %in% resultsNames(dds_tmp)) {
      res_shr <- lfcShrink(dds_tmp, coef = coef_name2, type = "apeglm")
      res_shr_df <- as.data.frame(res_shr); res_shr_df$gene <- rownames(res_shr_df)
      write.csv(res_shr_df, file.path(OUT_DIR, "tables", paste0("DE_", tag, "_LFCshrunk_apeglm.csv")), row.names = FALSE)
      ma_plot(res_shr, "apeglm_relevel")
    } else if (requireNamespace("ashr", quietly = TRUE)) {
      # Optional fallback: ashr supports contrast=
      res_shr <- lfcShrink(dds, contrast = c(var, levelA, levelB), type = "ashr")
      res_shr_df <- as.data.frame(res_shr); res_shr_df$gene <- rownames(res_shr_df)
      write.csv(res_shr_df, file.path(OUT_DIR, "tables", paste0("DE_", tag, "_LFCshrunk_ashr.csv")), row.names = FALSE)
      ma_plot(res_shr, "ashr")
    } else {
      message("Shrinkage skipped for ", tag, " (no direct coef; install 'ashr' for contrast shrinkage or accept unshrunk LFCs).")
    }
  }

  # Console summary + MA of unshrunk results
  sig <- subset(res_df, !is.na(padj) & padj < 0.05)
  sig_lfc1 <- subset(sig, abs(log2FoldChange) >= 1)
  message(sprintf("[%s] padj<0.05: %d; padj<0.05 & |LFC|>=1: %d", tag, nrow(sig), nrow(sig_lfc1)))

  png(file.path(OUT_DIR, "figs", paste0("de_02_MA_", tag, "_unshrunk.png")), width=1200, height=900, res=150)
  plotMA(res, ylim = c(-4,4), main = paste("MA (unshrunk):", tag))
  dev.off()

  invisible(list(res = res_df, res_shr = res_shr_df))
}

# 3.3 Run the contrasts you care about
out_MG_vs_E  <- run_contrast("condition_group", "SPACEFLIGHT_MICROGRAVITY", "EARTH",           "MG_vs_EARTH")
out_1G_vs_E  <- run_contrast("condition_group", "SPACEFLIGHT_1G",           "EARTH",           "1G_vs_EARTH")
out_MG_vs_1G <- run_contrast("condition_group", "SPACEFLIGHT_MICROGRAVITY", "SPACEFLIGHT_1G",  "MG_vs_1G")

# 3.6 Quick overall counts
cat("\n=== Overall DE summaries (padj < 0.05) ===\n")
for (nm in c("MG_vs_EARTH","1G_vs_EARTH","MG_vs_1G")) {
  f <- file.path(OUT_DIR, "tables", paste0("DE_", nm, "_unshrunk.csv"))
  if (!file.exists(f)) { 
    cat(sprintf("%-12s  (no file)\n", nm))
    next
  }
  df <- read.csv(f, check.names = FALSE)
  up   <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange > 0)
  down <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange < 0)
  cat(sprintf("%-12s  up: %5d   down: %5d   total: %5d\n", nm, up, down, up+down))
}

# 3.7 Tiny text summary
writeLines(c(
  paste("Design:", deparse(design_formula)),
  sprintf("MG vs EARTH  : %d sig (padj<0.05)", 
          sum(!is.na(out_MG_vs_E$res$padj)  & out_MG_vs_E$res$padj  < 0.05)),
  sprintf("1G vs EARTH  : %d sig (padj<0.05)", 
          sum(!is.na(out_1G_vs_E$res$padj)  & out_1G_vs_E$res$padj  < 0.05)),
  sprintf("MG vs 1G     : %d sig (padj<0.05)", 
          sum(!is.na(out_MG_vs_1G$res$padj) & out_MG_vs_1G$res$padj < 0.05))
), con = file.path(OUT_DIR, "tables", "DE_step3_summary.txt"))
```
# Visualization
```
# Load packages (install once if needed)
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
}
if (!requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Dm.eg.db", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
})

# Safety / dirs 
stopifnot(exists("vsd"), exists("meta"))       # from Step 3
stopifnot(exists("out_MG_vs_E"), exists("out_1G_vs_E"), exists("out_MG_vs_1G"))
if (!exists("OUT_DIR")) OUT_DIR <- "RESULTS_OSD514"
dir.create(file.path(OUT_DIR, "figs"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

# FBgn -> official symbol map (fallback to ID if missing)
fbgn_keys <- unique(rownames(vsd))
map_tbl <- AnnotationDbi::select(org.Dm.eg.db,
                                 keys = fbgn_keys,
                                 columns = "SYMBOL",
                                 keytype = "FLYBASE")
fbgn_sym_map <- setNames(map_tbl$SYMBOL, map_tbl$FLYBASE)
toSym <- function(x) {
  y <- unname(fbgn_sym_map[x])
  y[is.na(y) | y == ""] <- x[is.na(y) | y == ""]
  y
}

# Short column labels (E/1G/MG + sex initial + rep + lane)
if (!exists("short_col_lab")) {
  cond_code <- c(EARTH="E", SPACEFLIGHT_1G="1G", SPACEFLIGHT_MICROGRAVITY="MG")
  lab <- paste(cond_code[as.character(meta$condition_group)],
               substr(as.character(meta$sex), 1, 1),
               meta$replicate,
               meta$lane,
               sep = "_")
  names(lab) <- rownames(meta)
  short_col_lab <- lab[colnames(vsd)]
}

# Volcano plot helper

# Colors: Down=blue, NS=grey, Up=red; uses padj (FDR) 0.05 and |LFC|>=1 bands
save_volcano <- function(res_unshr, res_shr, tag, thr_p=0.05, thr_fc=1) {
  vdf <- as.data.frame(res_unshr)
  vdf$gene_id    <- rownames(vdf)
  vdf$log2FC_shr <- if (!is.null(res_shr)) res_shr[rownames(vdf), "log2FoldChange"] else vdf$log2FoldChange
  vdf$gene_symbol <- toSym(vdf$gene_id)

  vdf$group <- "NS"
  vdf$group[vdf$padj < thr_p & vdf$log2FC_shr >=  thr_fc] <- "Up (sig)"
  vdf$group[vdf$padj < thr_p & vdf$log2FC_shr <= -thr_fc] <- "Down (sig)"
  vdf$group <- factor(vdf$group, levels=c("Down (sig)","NS","Up (sig)"))

  # label top 15 by FDR (you can swap for a pathway set if you prefer)
  lab_set <- head(vdf$gene_id[order(vdf$padj)], 15)

  out <- file.path(OUT_DIR, "figs", paste0("volcano_", tag, ".png"))
  png(out, width=1800, height=1500, res=200, bg="white")
  print(
    ggplot(vdf, aes(x=log2FC_shr, y=-log10(padj), color=group)) +
      geom_point(alpha=0.7, size=1.6, na.rm=TRUE) +
      geom_vline(xintercept=c(-thr_fc, thr_fc), linetype="dashed") +
      geom_hline(yintercept=-log10(thr_p), linetype="dashed") +
      ggrepel::geom_text_repel(
        data = subset(vdf, gene_id %in% lab_set),
        aes(label=gene_symbol),
        size=3, max.overlaps=40, box.padding=0.4, min.segment.length=0
      ) +
      scale_color_manual(values=c("Down (sig)"="#2C7BB6",
                                  "NS"="grey60",
                                  "Up (sig)"="#D7191C")) +
      labs(x="Shrunken log2 fold change",
           y="-log10(FDR)",
           title=paste0("Volcano: ", gsub("_"," ", tag)),
           subtitle="Dashed lines: |log2FC|=1 and FDR=0.05") +
      theme_classic() + theme(legend.position="top")
  )
  dev.off()
}

# Heatmap helper
# Subsets to the two groups being compared; orders columns condition → sex → replicate.
save_heatmap <- function(tag, res_unshr, groups, topN=40, cap_z=2.5) {
  mat <- assay(vsd)

  # pick top genes by FDR then |LFC|
  rdf <- as.data.frame(res_unshr)
  rdf$gene_id <- rownames(rdf)
  rdf <- rdf[order(rdf$padj, -abs(rdf$log2FoldChange)), ]
  genes_use <- head(stats::na.omit(rdf$gene_id), topN)
  genes_use <- intersect(genes_use, rownames(mat))
  if (!length(genes_use)) stop("No genes available for heatmap: ", tag)

  # subset columns to requested groups
  cols_use <- rownames(meta)[meta$condition_group %in% groups]
  if (!length(cols_use)) stop("No columns match groups for heatmap: ", paste(groups, collapse=", "))
  plot_mat <- mat[genes_use, cols_use, drop=FALSE]

  # row-Z scale and cap
  z <- t(scale(t(plot_mat))); z[is.na(z)] <- 0
  z[z >  cap_z] <-  cap_z
  z[z < -cap_z] <- -cap_z
  rownames(z) <- toSym(rownames(z))  # FBgn -> symbol

  # annotations (only for selected columns)
  if (!"batch" %in% colnames(meta)) meta$batch <- factor(paste(meta$flowcell, meta$lane, sep="_"))
  ann <- meta[cols_use, c("condition_group","batch","sex"), drop=FALSE]
  # ensure replicate exists for ordering
  ann$replicate <- if ("replicate" %in% colnames(meta)) meta[cols_use, "replicate", drop=TRUE] else seq_len(nrow(ann))

  # explicit order: condition → sex → replicate
  cond_order_all <- c("EARTH","SPACEFLIGHT_1G","SPACEFLIGHT_MICROGRAVITY")
  cond_order <- cond_order_all[cond_order_all %in% groups]
  ord <- order(factor(ann$condition_group, levels=cond_order),
               as.character(ann$sex),
               as.numeric(ann$replicate))
  z   <- z[, ord, drop=FALSE]
  ann <- ann[ord, , drop=FALSE]
  lab_col <- short_col_lab[colnames(z)]

  out <- file.path(OUT_DIR, "figs", paste0("heatmap_", tag, ".png"))
  png(out, width=2200, height=1500, res=200, bg="white")
  pheatmap::pheatmap(z,
    scale="none",
    cluster_rows=TRUE,
    cluster_cols=FALSE,             # keep explicit order
    annotation_col=ann,
    labels_col=lab_col,
    show_rownames=TRUE,
    show_colnames=TRUE,
    fontsize_row=7,
    fontsize_col=9,
    border_color=NA,
    main=paste0("Top ", nrow(z), " DE genes (", gsub("_"," ", tag), "); row-Z cap ±", cap_z)
  )
  dev.off()
}

# Make plots
# Volcanoes (padj < 0.05; |LFC| ≥ 1 for color bands)
save_volcano(out_MG_vs_E$res,  out_MG_vs_E$res_shr,  "MG_vs_EARTH")
save_volcano(out_1G_vs_E$res,  out_1G_vs_E$res_shr,  "1G_vs_EARTH")
save_volcano(out_MG_vs_1G$res, out_MG_vs_1G$res_shr, "MG_vs_1G")

# Heatmaps (each contrast shows only its two groups)
save_heatmap("MG_vs_EARTH", out_MG_vs_E$res,
             groups=c("SPACEFLIGHT_MICROGRAVITY","EARTH"),
             topN=40, cap_z=2.5)

save_heatmap("1G_vs_EARTH", out_1G_vs_E$res,
             groups=c("SPACEFLIGHT_1G","EARTH"),
             topN=40, cap_z=2.5)

save_heatmap("MG_vs_1G", out_MG_vs_1G$res,
             groups=c("SPACEFLIGHT_MICROGRAVITY","SPACEFLIGHT_1G"),
             topN=40, cap_z=2.5)

message("Saved PNGs to: ", file.path(OUT_DIR, "figs"))
```
# Mito-Specific Analysis (GO ORA)

After analyzing and visualizing DESeq2 results, we decided to zero in on mitochondrially relevant genes and where they land on the spectrum of differentially expressed genes by performing GO Overrepresentation Analysis (ORA) on pathways of interest and related terms.

```
# Working dir (edit if needed)
setwd("~/Desktop/ADBR Mito/OSD-514")

# Load Packages
pkgs_cran <- c("ggplot2","ggrepel","pheatmap","dplyr")
for (p in pkgs_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("DESeq2","AnnotationDbi","org.Dm.eg.db","GO.db")
for (p in pkgs_bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(pheatmap); library(dplyr)
  library(DESeq2);  library(AnnotationDbi); library(org.Dm.eg.db); library(GO.db)
})

# Output dirs
if (!exists("OUT_DIR")) OUT_DIR <- "RESULTS_OSD514"
dir.create(file.path(OUT_DIR,"figs"),   recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR,"tables"), recursive=TRUE, showWarnings=FALSE)
tbl_dir <- file.path(OUT_DIR, "tables")

# Load DE results from prev step
load_de <- function(tag) {
  f_un  <- file.path(tbl_dir, paste0("DE_", tag, "_unshrunk.csv"))
  if (!file.exists(f_un)) stop("Missing unshrunk DE table: ", f_un)
  res_df <- read.csv(f_un, check.names=FALSE)
  if ("gene" %in% names(res_df)) rownames(res_df) <- res_df$gene

  # prefer apeglm files if present, otherwise ashr, else any file that starts with LFCshrunk_
  cand <- c("apeglm","apeglm_relevel","ashr")
  f_shr <- NULL
  for (m in cand) {
    f_try <- file.path(tbl_dir, paste0("DE_", tag, "_LFCshrunk_", m, ".csv"))
    if (file.exists(f_try)) { f_shr <- f_try; break }
  }
  if (is.null(f_shr)) {
    # fallback: pick any LFCshrunk_* file if exists
    f_any <- list.files(tbl_dir, pattern=paste0("^DE_", tag, "_LFCshrunk_.*\\.csv$"), full.names=TRUE)
    if (length(f_any)) f_shr <- f_any[1]
  }
  if (is.null(f_shr)) stop("No shrunken LFC table found for ", tag, ".")

  shr_df <- read.csv(f_shr, check.names=FALSE)
  if ("gene" %in% names(shr_df)) rownames(shr_df) <- shr_df$gene
  message("Loaded: ", basename(f_un), " + ", basename(f_shr))
  list(res=res_df, res_shr=shr_df)
}

out_MG_vs_E   <- load_de("MG_vs_EARTH")
out_1G_vs_E   <- load_de("1G_vs_EARTH")
out_MG_vs_1G  <- load_de("MG_vs_1G")

# Try to get vsd/meta (for heatmaps). Prefer DESeq2 RDS; else rebuild from Step 2 CSVs.
if (!exists("vsd") || !exists("meta")) {
  f_dds <- file.path(OUT_DIR,"tables","dds_step3_DESeq2.rds")
  if (file.exists(f_dds)) {
    dds  <- readRDS(f_dds)
    meta <- as.data.frame(colData(dds))
    vsd  <- vst(dds, blind=TRUE)
    message("Loaded vsd/meta from: ", f_dds)
  } else {
    f_counts <- file.path("Collapsed_Counts","tables","counts_DESeq_ready.csv")
    f_meta   <- file.path("Collapsed_Counts","tables","metadata_from_filenames.csv")
    if (file.exists(f_counts) && file.exists(f_meta)) {
      message("Rebuilding vsd/meta from Step-2 CSVs.")
      counts <- as.matrix(read.csv(f_counts, row.names=1, check.names=FALSE))
      meta   <- read.csv(f_meta,   row.names=1, check.names=FALSE)
      meta$condition_group <- factor(meta$condition_group,
        levels=c("EARTH","SPACEFLIGHT_1G","SPACEFLIGHT_MICROGRAVITY"))
      if (!is.factor(meta$sex)) meta$sex <- factor(meta$sex)
      # simple design for VST; DE already done
      dds <- DESeqDataSetFromMatrix(round(counts), colData=meta, design=~ sex + condition_group)
      dds <- estimateSizeFactors(dds)
      vsd <- vst(dds, blind=TRUE)
    } else {
      message("No vsd/meta available; heatmaps will be skipped (volcanoes still created).")
    }
  }
}

# FBgn → SYMBOL map (use union of DE IDs + vsd if present)
fbgn_keys <- unique(c(
  rownames(out_MG_vs_E$res), rownames(out_1G_vs_E$res), rownames(out_MG_vs_1G$res),
  if (exists("vsd")) rownames(vsd) else character(0)
))
map_tbl <- AnnotationDbi::select(org.Dm.eg.db, keys=fbgn_keys, columns="SYMBOL", keytype="FLYBASE")
fbgn_sym_map <- setNames(map_tbl$SYMBOL, map_tbl$FLYBASE)
toSym <- function(x) { y <- unname(fbgn_sym_map[x]); y[is.na(y) | y==""] <- x[is.na(y) | y==""]; y }

# Short sample labels for heatmaps
if (exists("meta") && exists("vsd") && !exists("short_col_lab")) {
  cond_code <- c(EARTH="E", SPACEFLIGHT_1G="1G", SPACEFLIGHT_MICROGRAVITY="MG")
  default_lane <- if ("lane" %in% colnames(meta)) meta$lane else ""
  lab <- paste(cond_code[as.character(meta$condition_group)],
               substr(as.character(meta$sex), 1, 1),
               if ("replicate" %in% colnames(meta)) meta$replicate else "",
               default_lane, sep="_")
  names(lab) <- rownames(meta)
  short_col_lab <- lab[colnames(vsd)]
}

# Curated mitochondrial GO sets (BP & CC + descendants)
mito_go_bp <- c(
  "GO:0006119",  # oxidative phosphorylation
  "GO:0022900",  # electron transport chain
  "GO:0006099",  # tricarboxylic acid (TCA) cycle
  "GO:0006635",  # fatty acid beta-oxidation
  "GO:0000422"   # mitophagy
)
mito_go_cc <- c(
  "GO:0005739",  # mitochondrion
  "GO:0005743",  # mitochondrial inner membrane
  "GO:0005747",  # respiratory chain complex I
  "GO:0005753",  # ATP synthase complex
  "GO:0005759"   # mitochondrial matrix
)

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

# Volcano (label curated mito genes; colors: Down=blue, NS=grey, Up=red)
save_volcano_mito <- function(res_unshr, res_shr, tag, mito_fbgn, thr_p=0.05, thr_fc=1) {
  vdf <- as.data.frame(res_unshr)
  vdf$gene_id     <- rownames(vdf)
  vdf$log2FC_shr  <- res_shr[rownames(vdf), "log2FoldChange"]
  vdf$gene_symbol <- toSym(vdf$gene_id)

  vdf$group <- "NS"
  vdf$group[vdf$padj < thr_p & vdf$log2FC_shr >=  thr_fc] <- "Up (sig)"
  vdf$group[vdf$padj < thr_p & vdf$log2FC_shr <= -thr_fc] <- "Down (sig)"
  vdf$group <- factor(vdf$group, levels=c("Down (sig)","NS","Up (sig)"))

  lab_set <- intersect(mito_fbgn, vdf$gene_id)
  if (!length(lab_set)) lab_set <- head(vdf[order(vdf$padj), "gene_id"], 12)

  out <- file.path(OUT_DIR,"figs", paste0("volcano_mito_", tag, ".png"))
  png(out, width=1800, height=1500, res=200, bg="white")
  print(
    ggplot(vdf, aes(x=log2FC_shr, y=-log10(padj), color=group)) +
      geom_point(alpha=0.7, size=1.6, na.rm=TRUE) +
      # halo the labeled mito genes
      geom_point(data=subset(vdf, gene_id %in% lab_set), shape=21, stroke=0.7, size=2.8, color="black") +
      geom_vline(xintercept=c(-thr_fc, thr_fc), linetype="dashed") +
      geom_hline(yintercept=-log10(thr_p), linetype="dashed") +
      ggrepel::geom_text_repel(
        data=subset(vdf, gene_id %in% lab_set),
        aes(label=gene_symbol),
        size=3, max.overlaps=40, box.padding=0.4, min.segment.length=0
      ) +
      scale_color_manual(values=c("Down (sig)"="#2C7BB6","NS"="grey60","Up (sig)"="#D7191C")) +
      labs(x="Shrunken log2 fold change", y="-log10(FDR)",
           title=paste0("Volcano (mitochondrial focus): ", gsub("_"," ", tag)),
           subtitle="Dashed: |log2FC|=1 and FDR=0.05. Labels = curated mito genes (GO).") +
      theme_classic() + theme(legend.position="top")
  )
  dev.off()
}

# Heatmap of curated mito genes (only the two groups being compared)
save_mito_heatmap <- function(tag, res_unshr, groups, mito_fbgn, max_rows=50, cap_z=2.5) {
  if (!exists("vsd") || !exists("meta")) { message("No vsd/meta; skipping heatmap for ", tag); return(invisible(NULL)) }
  mat <- assay(vsd)

  # choose the columns for the two groups
  cols_use <- rownames(meta)[meta$condition_group %in% groups]
  if (!length(cols_use)) { message("No columns match groups for ", tag); return(invisible(NULL)) }

  ann_cols <- intersect(c("condition_group","sex","replicate","flowcell","lane","batch"),
                        colnames(meta))
  ann <- meta[cols_use, ann_cols, drop=FALSE]

  # explicit order: condition -> sex -> replicate
  cond_levels <- c("EARTH","SPACEFLIGHT_1G","SPACEFLIGHT_MICROGRAVITY")
  ord_keys <- list(factor(ann$condition_group, levels=cond_levels[cond_levels %in% groups]))
  if ("sex" %in% colnames(ann))       ord_keys[[length(ord_keys)+1]] <- ann$sex
  if ("replicate" %in% colnames(ann)) ord_keys[[length(ord_keys)+1]] <- ann$replicate
  ord <- do.call(order, ord_keys)
  ann <- ann[ord, , drop=FALSE]
  cols_use <- rownames(ann)

  # pick curated mito genes present; rank by padj then |LFC|
  rdf <- as.data.frame(res_unshr); rdf$gene_id <- rownames(rdf)
  genes_in <- intersect(mito_fbgn, rownames(mat))
  if (!length(genes_in)) { message("Curated mito set has no overlap for ", tag); return(invisible(NULL)) }
  # guard if some genes missing from res_unshr (should be rare)
  genes_rankable <- intersect(genes_in, rownames(rdf))
  ord_idx <- order(rdf[genes_rankable, "padj"], -abs(rdf[genes_rankable, "log2FoldChange"]))
  genes_plot <- genes_rankable[ord_idx][seq_len(min(length(genes_rankable), max_rows))]

  z <- t(scale(t(mat[genes_plot, cols_use, drop=FALSE])))
  z[is.na(z)] <- 0
  z[z >  cap_z] <-  cap_z
  z[z < -cap_z] <- -cap_z
  rownames(z) <- toSym(rownames(z))

  if (!exists("short_col_lab")) short_col_lab <- setNames(colnames(z), colnames(z))
  lab_col <- short_col_lab[colnames(z)]

  # batch token if not present
  if (!"batch" %in% colnames(ann) && all(c("flowcell","lane") %in% colnames(ann))) {
    ann$batch <- factor(paste(ann$flowcell, ann$lane, sep="_"))
  }

  out <- file.path(OUT_DIR,"figs", paste0("heatmap_mito_", tag, ".png"))
  png(out, width=2200, height=1500, res=200, bg="white")
  pheatmap::pheatmap(z,
    scale="none",
    cluster_rows=TRUE,
    cluster_cols=FALSE,  # keep explicit order (no gaps)
    annotation_col = ann[, intersect(c("condition_group","sex","batch"), colnames(ann)), drop=FALSE],
    labels_col = lab_col,
    show_rownames=TRUE, show_colnames=TRUE,
    fontsize_row=8, fontsize_col=9,
    border_color=NA,
    main=paste0("Curated mitochondrial genes (", gsub("_"," ", tag), "); row-Z cap ±", cap_z,
                " (", nrow(z), " genes)")
  )
  dev.off()
}

# Make the plots
save_volcano_mito(out_MG_vs_E$res,  out_MG_vs_E$res_shr,  "MG_vs_EARTH", mito_fbgn_curated)
save_volcano_mito(out_1G_vs_E$res,  out_1G_vs_E$res_shr,  "1G_vs_EARTH", mito_fbgn_curated)
save_volcano_mito(out_MG_vs_1G$res, out_MG_vs_1G$res_shr, "MG_vs_1G",   mito_fbgn_curated)

save_mito_heatmap("MG_vs_EARTH",
  res_unshr = out_MG_vs_E$res,
  groups    = c("SPACEFLIGHT_MICROGRAVITY","EARTH"),
  mito_fbgn = mito_fbgn_curated,
  max_rows  = 50, cap_z=2.5)

save_mito_heatmap("1G_vs_EARTH",
  res_unshr = out_1G_vs_E$res,
  groups    = c("SPACEFLIGHT_1G","EARTH"),
  mito_fbgn = mito_fbgn_curated,
  max_rows  = 50, cap_z=2.5)

save_mito_heatmap("MG_vs_1G",
  res_unshr = out_MG_vs_1G$res,
  groups    = c("SPACEFLIGHT_MICROGRAVITY","SPACEFLIGHT_1G"),
  mito_fbgn = mito_fbgn_curated,
  max_rows  = 50, cap_z=2.5)

message("Saved focused PNGs to: ", file.path(OUT_DIR, "figs"))
```
# GSEA (in R)

```
# Working directory
setwd("~/Desktop/ADBR Mito/OSD-514")

# Packages
pkgs_cran <- c("data.table","ggplot2")
for (p in pkgs_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("fgsea","AnnotationDbi","org.Dm.eg.db","GO.db")
for (p in pkgs_bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(fgsea);      library(AnnotationDbi)
  library(org.Dm.eg.db); library(GO.db)
})

# Output dirs
GSEA_DIR <- "GSEA"
dir.create(GSEA_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(GSEA_DIR,"tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(GSEA_DIR,"figs"),   recursive=TRUE, showWarnings=FALSE)

OUT_DIR <- "RESULTS_OSD514"
TBL_DIR <- file.path(OUT_DIR,"tables")

# Helpers
load_de <- function(tag) {
  f <- file.path(TBL_DIR, paste0("DE_", tag, "_unshrunk.csv"))
  if (!file.exists(f)) stop("Missing DE table: ", f)
  df <- read.csv(f, check.names=FALSE)
  if ("gene" %in% names(df)) rownames(df) <- df$gene
  df
}

# Build ranks: prefer Wald stat; else signed -log10(p)
build_ranks <- function(res_df) {
  ids <- rownames(res_df)
  if ("stat" %in% names(res_df) && any(is.finite(res_df$stat))) {
    score <- res_df$stat
  } else {
    pv    <- if ("pvalue" %in% names(res_df)) res_df$pvalue else res_df$padj
    pv[!is.finite(pv) | is.na(pv)] <- 1
    lfc   <- if ("log2FoldChange" %in% names(res_df)) res_df$log2FoldChange else 0
    score <- sign(lfc) * (-log10(pmax(pv, 1e-300)))
  }
  keep <- is.finite(score) & !is.na(ids) & nzchar(ids)
  ranks <- tapply(score[keep], ids[keep], mean)  # aggregate duplicates
  sort(ranks, decreasing=TRUE)
}

# Build GO BP sets limited to the current ranked universe (robust)
build_go_bp_sets <- function(ranks, min_size=10, max_size=500) {
  genes <- names(ranks)
  m <- AnnotationDbi::select(org.Dm.eg.db,
                             keys=genes, keytype="FLYBASE",
                             columns=c("GOALL","ONTOLOGYALL"))
  m <- unique(as.data.frame(m, stringsAsFactors=FALSE))
  m <- m[!is.na(m$GOALL) & !is.na(m$ONTOLOGYALL) & m$ONTOLOGYALL=="BP", , drop=FALSE]
  pathways <- split(m$FLYBASE, m$GOALL, drop=TRUE)
  pathways <- lapply(pathways, function(x) unique(x[!is.na(x) & nzchar(x)]))
  lens <- vapply(pathways, length, integer(1))
  pathways[lens >= min_size & lens <= max_size]
}

# Map GOID -> TERM
term_map_for <- function(go_ids) {
  if (!length(go_ids)) return(setNames(character(0), character(0)))
  tbl <- AnnotationDbi::select(GO.db, keys=go_ids, keytype="GOID", columns="TERM")
  tbl <- unique(as.data.frame(tbl, stringsAsFactors=FALSE))
  setNames(tbl$TERM, tbl$GOID)
}

# Null-coalesce
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

# Pick mito-relevant pathways; fallback to top hits if none
pick_mito_paths <- function(fg_dt, padj_cut=0.25, top_n=8) {
  if (!nrow(fg_dt)) return(character(0))
  lab <- if ("term" %in% names(fg_dt)) fg_dt$term else fg_dt$pathway
  mito_pat <- paste(
    "mitochond", "oxidative phosph", "electron transport", "respiratory chain",
    "ATP synthase", "tricarboxylic", "TCA", "beta-oxid", "mitophagy", "fission", "fusion",
    sep="|"
  )
  mito <- fg_dt[grepl(mito_pat, lab, ignore.case=TRUE) & is.finite(fg_dt$padj), , drop=FALSE]
  mito <- mito[order(mito$padj, -abs(mito$NES)), , drop=FALSE]
  sel <- mito$pathway[mito$padj <= padj_cut]
  if (!length(sel)) sel <- head(mito$pathway, min(top_n, nrow(mito)))
  unique(na.omit(sel))
}

# Save enrichment curve for one pathway
save_enrichment_curve <- function(tag, pathway_id, ranks, pathways, TERM_MAP) {
  if (!pathway_id %in% names(pathways)) return(invisible(NULL))
  term_label <- TERM_MAP[[pathway_id]] %||% pathway_id
  gp <- fgsea::plotEnrichment(pathways[[pathway_id]], ranks) +
        ggplot2::ggtitle(paste0(term_label, " (", tag, ")"))
  fn <- file.path(GSEA_DIR,"figs",
                  paste0("fgsea_curve_",
                         gsub("[^A-Za-z0-9]+","_", term_label), "_", tag, ".png"))
  ggsave(fn, gp, width=7, height=5, dpi=300, bg="white")
  message("Saved: ", fn)
}

# Save compact GSEA table figure (strip with NES/FDR)
save_gsea_table_plot <- function(tag, sel_paths, ranks, pathways, fg_dt, TERM_MAP) {
  sel_paths <- intersect(sel_paths, names(pathways))
  if (!length(sel_paths)) { message("No selected pathways to plot for ", tag); return(invisible(NULL)) }
  fg_sel <- fg_dt[match(sel_paths, fg_dt$pathway), , drop=FALSE]
  fn <- file.path(GSEA_DIR, "figs", paste0("fgsea_table_", tag, ".png"))
  png(fn, width=1600, height=900, res=150, bg="white")
  fgsea::plotGseaTable(pathways[sel_paths], ranks, fg_sel)
  dev.off()
  message("Saved: ", fn)
}

# Run everything for one contrast
run_fgsea_and_plots <- function(tag, nperm=10000, padj_cut_table=0.25, topN_bar=15) {
  message("\n=== ", tag, " ===")
  res_df <- load_de(tag)
  ranks  <- build_ranks(res_df)
  if (length(ranks) < 50) stop("Too few ranked genes for ", tag)

  pathways <- build_go_bp_sets(ranks, min_size=10, max_size=500)
  if (!length(pathways)) stop("No GO BP pathways built for ", tag)
  TERM_MAP <- term_map_for(names(pathways))

  set.seed(42)
  fg <- suppressWarnings(
    fgsea(pathways=pathways, stats=ranks, minSize=10, maxSize=500, nperm=nperm)
  )
  fg_dt <- as.data.table(fg)
  if (!nrow(fg_dt)) { message("fgsea returned no results for ", tag); return(invisible(NULL)) }
  setorder(fg_dt, padj, -NES)
  fg_dt[, term := TERM_MAP[pathway]]

  # Save full table
  out_csv <- file.path(GSEA_DIR, "tables", paste0("fgsea_GO_BP_", tag, ".csv"))
  fwrite(fg_dt[, .(pathway, term, size, NES, pval, padj,
                   leadingEdge = vapply(leadingEdge, \(v) paste(v, collapse=";"), character(1)))],
         out_csv)
  message("Saved: ", out_csv)

  # Quick top barplot
  topN <- head(fg_dt[is.finite(padj)], topN_bar)
  if (nrow(topN)) {
    topN[, label := ifelse(is.na(term) | term=="", pathway, term)]
    p <- ggplot(topN, aes(x=reorder(label, NES), y=NES, fill=-log10(padj))) +
      geom_col() + coord_flip() +
      labs(title=paste0("fgsea GO BP: ", gsub("_"," ", tag)),
           x=NULL, y="NES", fill="-log10(FDR)") +
      theme_minimal(base_size=12)
    out_png <- file.path(GSEA_DIR, "figs", paste0("fgsea_GO_BP_", tag, ".png"))
    ggsave(out_png, p, width=8, height=6, dpi=300, bg="white")
    message("Saved: ", out_png)
  }

  # Pick mito-relevant terms; if none, fall back so plots are still produced
  sel <- pick_mito_paths(fg_dt, padj_cut=padj_cut_table, top_n=8)
  if (!length(sel)) {
    message("No mito keywords found at the chosen FDR; falling back to top 6 terms for curves.")
    sel <- head(fg_dt$pathway, 6)
  }

  # Enrichment curves
  invisible(lapply(sel, save_enrichment_curve,
                   tag=tag, ranks=ranks, pathways=pathways, TERM_MAP=TERM_MAP))

  # Summary table figure
  save_gsea_table_plot(tag, sel, ranks, pathways, fg_dt, TERM_MAP)
}

# Run all three contrasts
for (tag in c("MG_vs_EARTH","1G_vs_EARTH","MG_vs_1G")) {
  run_fgsea_and_plots(tag)
}

message("\nAll done. See: ", normalizePath(GSEA_DIR))
```
