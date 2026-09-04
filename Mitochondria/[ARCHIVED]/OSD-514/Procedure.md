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
# GSEA
```
# Package install
pkgs_cran <- c(
  "data.table", "ggplot2", "dplyr", "ggrepel", "pheatmap",
  "enrichR", "stringr", "tidyr", "igraph", "ggraph"
)

pkgs_bioc <- c(
  "fgsea", "AnnotationDbi", "org.Dm.eg.db", "GO.db"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  library(GO.db)
  library(enrichR)
  library(stringr)
  library(tidyr)
  library(igraph)
  library(ggraph)
})

# Parameters
OUT_DIR   <- "RESULTS_OSD514"
TBL_DIR   <- file.path(OUT_DIR, "tables")
GSEA_DIR  <- "GSEA"

PADJ_CUT  <- 0.05
LFC_CUT   <- 1
FGSEA_NPERM <- 10000
TOP_BP_TERMS <- 15
TOP_MITO_CNET_TERMS <- 12

CONTRASTS <- c("MG_vs_EARTH", "1G_vs_EARTH", "MG_vs_1G")

dir.create(GSEA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(GSEA_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(GSEA_DIR, "figs"),   recursive = TRUE, showWarnings = FALSE)

# Helpers
`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0) a else b
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

wrap_terms <- function(x, width = 40) {
  stringr::str_wrap(x, width = width)
}

# Load DE results
load_de <- function(tag) {
  f <- file.path(TBL_DIR, paste0("DE_", tag, "_unshrunk.csv"))
  if (!file.exists(f)) stop("Missing DE table: ", f)

  df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)

  if ("gene" %in% colnames(df)) {
    rownames(df) <- df$gene
  }

  if (is.null(rownames(df)) || any(rownames(df) == "")) {
    stop("No usable gene identifiers found for ", tag, ". Expected 'gene' column or rownames.")
  }

  df
}

# Build fgsea ranks
# Prefer Wald stat; fallback to signed -log10(p)
build_ranks <- function(res_df) {
  ids <- rownames(res_df)

  if ("stat" %in% names(res_df) && any(is.finite(res_df$stat), na.rm = TRUE)) {
    score <- res_df$stat
  } else {
    pv <- if ("pvalue" %in% names(res_df)) res_df$pvalue else res_df$padj
    pv[!is.finite(pv) | is.na(pv)] <- 1

    lfc <- if ("log2FoldChange" %in% names(res_df)) res_df$log2FoldChange else 0
    lfc[!is.finite(lfc) | is.na(lfc)] <- 0

    score <- sign(lfc) * (-log10(pmax(pv, 1e-300)))
  }

  keep <- is.finite(score) & !is.na(ids) & nzchar(ids)
  ranks <- tapply(score[keep], ids[keep], mean)
  sort(ranks, decreasing = TRUE)
}

# Significant genes for ORA
get_sig_gene_ids <- function(res_df, padj_cut = 0.05, lfc_cut = 1) {
  req_cols <- c("padj", "log2FoldChange")
  miss <- setdiff(req_cols, colnames(res_df))
  if (length(miss)) stop("Missing required columns for ORA: ", paste(miss, collapse = ", "))

  df <- res_df[
    !is.na(res_df$padj) &
      is.finite(res_df$padj) &
      res_df$padj < padj_cut &
      !is.na(res_df$log2FoldChange) &
      is.finite(res_df$log2FoldChange) &
      abs(res_df$log2FoldChange) >= lfc_cut,
    ,
    drop = FALSE
  ]

  unique(rownames(df))
}

# FlyBase -> SYMBOL
fbgn_to_symbol <- function(ids) {
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (!length(ids)) return(character(0))

  tbl <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = ids,
    columns = "SYMBOL",
    keytype = "FLYBASE"
  )

  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  tbl <- tbl[!duplicated(tbl$FLYBASE), , drop = FALSE]

  map <- setNames(tbl$SYMBOL, tbl$FLYBASE)
  syms <- unname(map[ids])
  syms[is.na(syms) | syms == ""] <- ids[is.na(syms) | syms == ""]
  unique(syms)
}

# GO BP gene sets for fgsea
build_go_bp_sets <- function(ranks, min_size = 10, max_size = 500) {
  genes <- names(ranks)

  m <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = genes,
    keytype = "FLYBASE",
    columns = c("GOALL", "ONTOLOGYALL")
  )

  m <- unique(as.data.frame(m, stringsAsFactors = FALSE))
  m <- m[!is.na(m$GOALL) & !is.na(m$ONTOLOGYALL) & m$ONTOLOGYALL == "BP", , drop = FALSE]

  pathways <- split(m$FLYBASE, m$GOALL, drop = TRUE)
  pathways <- lapply(pathways, function(x) unique(x[!is.na(x) & nzchar(x)]))

  lens <- vapply(pathways, length, integer(1))
  pathways[lens >= min_size & lens <= max_size]
}

# GO term map
get_go_terms <- function(go_ids) {
  if (!length(go_ids)) return(setNames(character(0), character(0)))

  tbl <- AnnotationDbi::select(
    GO.db,
    keys = go_ids,
    keytype = "GOID",
    columns = "TERM"
  )

  tbl <- unique(as.data.frame(tbl, stringsAsFactors = FALSE))
  setNames(tbl$TERM, tbl$GOID)
}

# Mito keyword matcher
is_mito_term <- function(x) {
  grepl(
    paste(
      "mitochond",
      "oxidative phosph",
      "electron transport",
      "respiratory chain",
      "ATP synthase",
      "tricarboxylic",
      "TCA",
      "beta-oxid",
      "mitophagy",
      "fission",
      "fusion",
      sep = "|"
    ),
    x,
    ignore.case = TRUE
  )
}

# Save fgsea barplot
save_fgsea_barplot <- function(fg_dt, tag, top_n = 15) {

  top <- fg_dt[is.finite(padj) & !is.na(term), , drop = FALSE]
  top <- top[order(top$padj, -abs(top$NES)), , drop = FALSE]
  top <- head(top, top_n)

  if (!nrow(top)) return(invisible(NULL))

  top$label <- ifelse(is.na(top$term) | top$term == "", top$pathway, top$term)

  # Wrap long GO terms
  top$label_wrapped <- wrap_terms(top$label, width = 45)

  # Dynamic height based on number of terms
  plot_height <- max(6, 0.45 * nrow(top))

  p <- ggplot(top, aes(x = reorder(label_wrapped, NES), y = NES, fill = -log10(padj))) +
    geom_col(width = 0.8) +
    coord_flip(clip = "off") +
    scale_fill_viridis_c(option = "plasma", direction = -1) +
    labs(
      title = paste0("fgsea GO BP: ", gsub("_", " ", tag)),
      x = NULL,
      y = "Normalized Enrichment Score (NES)",
      fill = "-log10(FDR)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10, lineheight = 0.95),
      plot.margin = margin(10, 30, 10, 10),  # extra right margin for labels
      legend.position = "right"
    )

  out_png <- file.path(GSEA_DIR, "figs", paste0("fgsea_GO_BP_", tag, ".png"))

  ggsave(
    out_png,
    p,
    width = 10,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )
}

# Save fgsea enrichment curves
save_fgsea_curves <- function(fg_dt, pathways, ranks, term_map, tag, n_terms = 6) {
  mito <- fg_dt[!is.na(term) & is_mito_term(term) & is.finite(padj), , drop = FALSE]
  mito <- mito[order(mito$padj, -abs(mito$NES)), , drop = FALSE]

  sel <- head(mito$pathway, n_terms)

  if (!length(sel)) {
    sel <- head(fg_dt$pathway[order(fg_dt$padj, -abs(fg_dt$NES))], n_terms)
  }

  sel <- unique(sel[sel %in% names(pathways)])

  for (pw in sel) {
    term_label <- term_map[[pw]] %||% pw

    gp <- fgsea::plotEnrichment(pathways[[pw]], ranks) +
      ggtitle(paste0(term_label, " (", tag, ")"))

    out_png <- file.path(
      GSEA_DIR, "figs",
      paste0("fgsea_curve_", safe_filename(term_label), "_", tag, ".png")
    )

    ggsave(out_png, gp, width = 7, height = 5, dpi = 300, bg = "white")
  }
}

# Run enrichR
run_enrichr <- function(tag, gene_symbols) {
  if (!length(gene_symbols)) {
    message("No significant genes for enrichR in ", tag)
    return(NULL)
  }

  setEnrichrSite("FlyEnrichr")
  fly_dbs <- listEnrichrDbs()
  dbs_to_use <- grep("^GO_", fly_dbs$libraryName, value = TRUE)

  if (!length(dbs_to_use)) {
    warning("No GO libraries found in FlyEnrichr.")
    return(NULL)
  }

  enrich <- enrichr(gene_symbols, dbs_to_use)

  # save all enrichR tables
  for (db in names(enrich)) {
    df <- enrich[[db]]
    if (is.null(df) || !nrow(df)) next

    out_csv <- file.path(GSEA_DIR, "tables", paste0(tag, "_", safe_filename(db), ".csv"))
    fwrite(as.data.table(df), out_csv)
  }

  enrich
}

# Save enrichR plot
save_enrichr_plot <- function(enrich_df, title, outfile) {

  if (is.null(enrich_df) || !nrow(enrich_df)) return(invisible(NULL))

  # Detect GO category
  is_mf <- grepl("Molecular", outfile, ignore.case = TRUE)

  # Adaptive settings
  wrap_width   <- if (is_mf) 40 else 30     # tighter wrapping for MF
  height_scale <- if (is_mf) 1.0 else 0.6  # more vertical space for MF

  df <- enrich_df %>%
    dplyr::filter(!is.na(Adjusted.P.value)) %>%
    dplyr::arrange(Adjusted.P.value) %>%
    head(15)

  if (!nrow(df)) return(invisible(NULL))

  # Wrap labels
  df$label_wrapped <- stringr::str_wrap(df$Term, width = wrap_width)

  # Dynamic height
  plot_height <- max(7, height_scale * nrow(df))

  # Plot
  p <- ggplot(
    df,
    aes(
      x = reorder(label_wrapped, -log10(Adjusted.P.value)),
      y = -log10(Adjusted.P.value),
      fill = Combined.Score
    )
  ) +
    geom_col(width = 0.7) +
    coord_flip(clip = "off") +

    scale_fill_viridis_c(option = "viridis") +

    labs(
      title = title,
      x = NULL,
      y = "-log10(FDR)",
      fill = "Combined Score"
    ) +

    theme_minimal(base_size = 12) +

    theme(
      axis.text.y = element_text(
        size = 10,
        lineheight = if (is_mf) 1.0 else 0.95,   # extra spacing for MF
        margin = margin(r = 10)
      ),

      legend.position = "right",          
      plot.margin = margin(10, 10, 25, 10)     # space for legend
    )

  ggsave(
    outfile,
    p,
    width = 10,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )
}

# Build mito gene list from SYMBOLs
# based on GO annotation containing "mitochond"
find_mito_genes_from_symbols <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) return(character(0))

  go <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = symbols,
    columns = c("GO", "ONTOLOGY", "SYMBOL"),
    keytype = "SYMBOL"
  )

  go <- as.data.frame(go, stringsAsFactors = FALSE)

  go_terms <- AnnotationDbi::select(
    GO.db,
    keys = unique(go$GO[!is.na(go$GO)]),
    columns = c("TERM"),
    keytype = "GOID"
  )

  go_terms <- as.data.frame(go_terms, stringsAsFactors = FALSE)

  go_annot <- dplyr::left_join(go, go_terms, by = c("GO" = "GOID"))

  mito_annotations <- go_annot %>%
    dplyr::filter(!is.na(TERM) & grepl("mitochond", TERM, ignore.case = TRUE))

  unique(mito_annotations$SYMBOL)
}

# Extract enriched BP terms containing mito genes
extract_mito_terms_from_enrichr <- function(enrich_bp, mito_genes) {
  if (is.null(enrich_bp) || !nrow(enrich_bp) || !length(mito_genes)) {
    return(data.frame())
  }

  enrich_bp %>%
    as_tibble() %>%
    mutate(
      gene_list = strsplit(Genes, ";"),
      mito_hits = sapply(gene_list, function(x) {
        paste(intersect(x, mito_genes), collapse = ";")
      })
    ) %>%
    filter(mito_hits != "")
}

# Save mito dot plot
save_mito_dotplot <- function(mito_terms, tag) {
  if (is.null(mito_terms) || !nrow(mito_terms)) return(invisible(NULL))

  p <- mito_terms %>%
    mutate(score = -log10(Adjusted.P.value)) %>%
    ggplot(aes(x = score, y = reorder(Term, score))) +
    geom_point(size = 4) +
    labs(
      x = "-log10(FDR)",
      y = "GO term",
      title = paste0("Enriched pathways containing mitochondrial genes: ", gsub("_", " ", tag))
    ) +
    theme_bw()

  out_png <- file.path(GSEA_DIR, "figs", paste0("mito_terms_dotplot_", tag, ".png"))
  ggsave(out_png, p, width = 10, height = 7, dpi = 300, bg = "white")
}

# Build mito cnet object
build_mito_cnet <- function(enrich_df, mito_genes, top_n_terms = 15, max_term_chars = 55) {
  term_gene_df <- enrich_df %>%
    as_tibble() %>%
    filter(!is.na(Term), !is.na(Adjusted.P.value), !is.na(Genes)) %>%
    mutate(
      gene_list = strsplit(Genes, ";"),
      mito_hits = lapply(gene_list, function(x) intersect(x, mito_genes)),
      n_mito_hits = lengths(mito_hits)
    ) %>%
    filter(n_mito_hits > 0) %>%
    arrange(Adjusted.P.value, desc(n_mito_hits)) %>%
    slice_head(n = top_n_terms) %>%
    transmute(
      term_full = Term,
      term = str_trunc(Term, max_term_chars),
      Adjusted.P.value,
      minus_log10_fdr = -log10(Adjusted.P.value),
      mito_hits
    ) %>%
    unnest_longer(mito_hits, values_to = "gene") %>%
    filter(!is.na(gene), gene != "")

  if (nrow(term_gene_df) == 0) {
    return(NULL)
  }

  edges <- term_gene_df %>%
    dplyr::select(from = term, to = gene, minus_log10_fdr, term_full)

  term_nodes <- term_gene_df %>%
    distinct(term, term_full, minus_log10_fdr) %>%
    transmute(
      name = term,
      node_type = "term",
      score = minus_log10_fdr,
      label = term
    )

  gene_nodes <- term_gene_df %>%
    distinct(gene) %>%
    transmute(
      name = gene,
      node_type = "gene",
      score = NA_real_,
      label = gene
    )

  nodes <- bind_rows(term_nodes, gene_nodes)

  list(edges = edges, nodes = nodes, term_gene_df = term_gene_df)
}

# Save mito cnet plot
save_mito_cnet_plot <- function(cnet_obj, tag) {
  if (is.null(cnet_obj)) return(invisible(NULL))

  g <- graph_from_data_frame(
    d = cnet_obj$edges %>% dplyr::select(from, to),
    vertices = cnet_obj$nodes,
    directed = FALSE
  )

  set.seed(123)

  p <- ggraph(g, layout = "fr") +
    geom_edge_link(
      color = "grey70",
      alpha = 0.5,
      linewidth = 0.7
    ) +
    geom_node_point(
      aes(
        size = ifelse(node_type == "term", score, 4),
        fill = score
      ),
      shape = 21,
      color = "black",
      stroke = 0.3
    ) +
    geom_node_text(
      aes(label = label),
      repel = TRUE,
      size = 4
    ) +
    scale_fill_gradient(
      low = "blue",
      high = "red",
      na.value = "gold",
      name = expression(-log[10]("FDR"))
    ) +
    scale_size_continuous(
      range = c(4, 10),
      guide = "none"
    ) +
    labs(
      title = paste0("Mitochondrial genes driving enriched pathways: ", gsub("_", " ", tag)),
      subtitle = "Blue → red indicates increasing pathway significance"
    ) +
    theme_void()

  out_png <- file.path(GSEA_DIR, "figs", paste0("mito_cnet_", tag, ".png"))
  ggsave(out_png, p, width = 10, height = 8, dpi = 300, bg = "white")
}

# Main per-contrast analysis
run_analysis <- function(tag) {
  message("\n==============================")
  message("Running: ", tag)
  message("==============================")

  # Load DE
  res <- load_de(tag)

  # fgsea
  ranks <- build_ranks(res)
  if (length(ranks) < 50) {
    warning("Too few ranked genes for fgsea in ", tag)
  } else {
    pathways <- build_go_bp_sets(ranks, min_size = 10, max_size = 500)
    term_map <- get_go_terms(names(pathways))

    if (length(pathways)) {
      set.seed(42)
      fg <- suppressWarnings(
        fgsea(
          pathways = pathways,
          stats = ranks,
          minSize = 10,
          maxSize = 500,
          nperm = FGSEA_NPERM
        )
      )

      fg_dt <- as.data.table(fg)
      if (nrow(fg_dt)) {
        fg_dt[, term := term_map[pathway]]

        out_csv <- file.path(GSEA_DIR, "tables", paste0("fgsea_GO_BP_", tag, ".csv"))
        fwrite(
          fg_dt[, .(
            pathway, term, size, NES, pval, padj,
            leadingEdge = vapply(leadingEdge, function(v) paste(v, collapse = ";"), character(1))
          )],
          out_csv
        )

        save_fgsea_barplot(fg_dt, tag, top_n = TOP_BP_TERMS)
        save_fgsea_curves(fg_dt, pathways, ranks, term_map, tag, n_terms = 6)
      }
    }
  }

  # significant genes for ORA
  sig_ids <- get_sig_gene_ids(res, padj_cut = PADJ_CUT, lfc_cut = LFC_CUT)
  sig_symbols <- fbgn_to_symbol(sig_ids)

  sig_table <- data.frame(
    flybase_id = sig_ids,
    symbol = fbgn_to_symbol(sig_ids),
    stringsAsFactors = FALSE
  )
  write.csv(
    sig_table,
    file.path(GSEA_DIR, "tables", paste0("sig_genes_", tag, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )

  # enrichR ORA
  enrich <- run_enrichr(tag, sig_symbols)

  if (!is.null(enrich)) {
    # BP
    if ("GO_Biological_Process_2018" %in% names(enrich)) {
      save_enrichr_plot(
        enrich[["GO_Biological_Process_2018"]],
        title = gsub("_", " ", tag),
        outfile = file.path(GSEA_DIR, "figs", paste0("GO_Biological_Process_", tag, ".png"))
      )
    }

    # MF
    if ("GO_Molecular_Function_2018" %in% names(enrich)) {
      save_enrichr_plot(
        enrich[["GO_Molecular_Function_2018"]],
        title = gsub("_", " ", tag),
        outfile = file.path(GSEA_DIR, "figs", paste0("GO_Molecular_Function_", tag, ".png"))
      )
    }

    # CC
    if ("GO_Cellular_Component_2018" %in% names(enrich)) {
      save_enrichr_plot(
        enrich[["GO_Cellular_Component_2018"]],
        title = gsub("_", " ", tag),
        outfile = file.path(GSEA_DIR, "figs", paste0("GO_Cellular_Component_", tag, ".png"))
      )
    }

    # mito gene extraction + mito term filtering from BP
    if ("GO_Biological_Process_2018" %in% names(enrich)) {
      mito_genes <- find_mito_genes_from_symbols(sig_symbols)

      write.csv(
        data.frame(gene = mito_genes, stringsAsFactors = FALSE),
        file.path(GSEA_DIR, "tables", paste0("mito_genes_", tag, ".csv")),
        row.names = FALSE,
        quote = FALSE
      )

      mito_terms <- extract_mito_terms_from_enrichr(
        enrich_bp = enrich[["GO_Biological_Process_2018"]],
        mito_genes = mito_genes
      )

      # NETWORK (robust)

    if ("GO_Biological_Process_2018" %in% names(enrich)) {

      enrich_bp <- enrich[["GO_Biological_Process_2018"]]

      if (!is.null(enrich_bp) && nrow(enrich_bp)) {

        # Step 1: try mito-based network
        mito_genes <- find_mito_genes_from_symbols(sig_symbols)

        message("Mito genes found: ", length(mito_genes))

        mito_terms <- extract_mito_terms_from_enrichr(enrich_bp, mito_genes)

        message("Mito-enriched terms: ", nrow(mito_terms))

        if (nrow(mito_terms) > 0) {

          message("Building mito-specific network")

          cnet_obj <- build_mito_cnet(
            enrich_df = enrich_bp,
            mito_genes = mito_genes,
            top_n_terms = 12
          )

        } else {

          # FALLBACK: use top enriched terms regardless of mito
          message("No mito terms found → using top enriched GO terms instead")

          top_terms <- enrich_bp %>%
            dplyr::filter(!is.na(Adjusted.P.value)) %>%
            dplyr::arrange(Adjusted.P.value) %>%
            head(12)

          cnet_obj <- build_mito_cnet(
            enrich_df = top_terms,
            mito_genes = unique(unlist(strsplit(top_terms$Genes, ";"))),
            top_n_terms = 12
          )
        }

        # Save plot if graph exists
        if (!is.null(cnet_obj)) {

          save_mito_cnet_plot(cnet_obj, tag)

        } else {
          message("Network could not be built for ", tag)
        }
      } else {
            message("No enriched GO BP terms containing mitochondrial genes for ", tag)
        }
      }
    }
  }
}   
# Run all contrasts
for (tag in CONTRASTS) {
  tryCatch(
    run_analysis(tag),
    error = function(e) {
      message("ERROR in ", tag, ": ", e$message)
    }
  )
}

message("\nDone. Outputs written to: ", normalizePath(GSEA_DIR))
```
