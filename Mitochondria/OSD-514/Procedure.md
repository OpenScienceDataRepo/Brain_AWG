# Extract Raw Counts

OSD-514 did not have a clean raw counts file to run in DESeq2. It had an "unnormalized counts" file but that had non-integer values, meaning it was a convenience matrix of RSEM “expected counts” (decimals).

The 24 files compiled are also not raw counts files. The *.genes.results files from RSEM contain expected counts (fractional). We used them because DESeq2 supports this workflow via tximport. We gave DESeq2 the expected counts plus the average transcript lengths that tximport extracts. DESeq2 then uses these as offsets (length-aware normalization).

Just rounding the “Unnormalized Counts” CSV throws away the information RSEM provides. It biases low counts, meaning small fractional differences can flip to 0 or 1 after rounding, distorting dispersion and p-values. It also breaks length-aware normalization since tximport+DESeq2 can properly correct for gene length and library size.

Because of this, the 24 raw count results were downloaded for each sample in each condition and collapsed into one matrix manually using the below code.

```
# Set Working Directory
setwd("~/Desktop/ADBR Mito/OSD-514")

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
