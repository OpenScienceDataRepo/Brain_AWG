setwd("/Volumes/Marians_SSD/ADBR_Mito/OSD-514/RNA_Seq/RawCounts")
# Set Working Directory

# Set the output folder
out_dir <- "/Volumes/Marians_SSD/ADBR_Mito/OSD-514/RNA_Seq/Collapsed_Counts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Find the 24 gene files
list.files()

files <- list.files(pattern="\\.genes\\.results$", full.names=TRUE)
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
write.csv(mat, file.path(out_dir, "OSD514_RSEM_expected_counts.csv"))

# Also write a rounded integer matrix that DESeq2 will accept directly
mat_int <- round(mat)
storage.mode(mat_int) <- "integer"
write.csv(mat_int, file.path(out_dir, "counts_DESeq_ready.csv"))

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
write.csv(meta, file.path(out_dir, "metadata_from_filenames.csv"))

cat("Wrote:\n",
    "- tables/OSD514_RSEM_expected_counts.csv\n",
    "- tables/counts_DESeq_ready.csv (rounded; plug into DESeq2)\n",
    "- tables/metadata_from_filenames.csv\n", sep = "")