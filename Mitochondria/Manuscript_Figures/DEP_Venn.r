# Mitochondrial DEP Venn Diagram (Proteomics)
# LFC threshold is 0.5 (|log2FC| >= 0.5), matching the OSD-514 volcano plot.
# The mito universe is matched both by gene symbol AND by UniProt AC to catch OSD-514 proteins that may be missing from org.Dm.eg.db's UNIPROT mapping.

# Packages
pkgs_cran <- c("ggplot2", "dplyr", "ggVennDiagram", "stringr")
pkgs_bioc <- c("AnnotationDbi", "org.Dm.eg.db", "GO.db")

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggVennDiagram)
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  library(GO.db)
})

# Paths
OSD514_FILES <- c(
  "Your_Working_Directory/Limma_SFug_vs_Earth_results_ms3fix.csv"
)

TAU_FILES <- c(
  "Your_Working_Directory/DE_TauR406W_vs_Control_Day20.csv"
)

OUT_DIR <- "Your_Working_Directory/CrossDataset_Venn_DEP"
dir.create(file.path(OUT_DIR, "figs"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

# Significance thresholds  (limma column names)
PADJ_CUT <- 0.05
LFC_CUT  <- 0.5    # |log2FC| >= 0.5, matching the OSD-514 volcano plot

PADJ_COL <- "adj.P.Val"
LFC_COL  <- "logFC"

# ID-parsing helpers
parse_uniprot_ac <- function(ids) {
  ifelse(
    grepl("\\|", ids),
    sub("^[a-z]+\\|([A-Z0-9]+)\\|.*$", "\\1", ids, perl = TRUE),
    ids
  )
}

# UniProt ACs -> data.frame: input_ac | symbol | mapped
uniprot_to_symbol_df <- function(acs) {
  acs <- unique(acs[!is.na(acs) & nzchar(acs)])
  if (!length(acs)) return(data.frame(input_ac=character(), symbol=character(), mapped=logical()))
  tbl <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys    = acs,
    keytype = "UNIPROT",
    columns = "SYMBOL"
  )
  tbl <- tbl[!duplicated(tbl$UNIPROT), ]
  map  <- setNames(tbl$SYMBOL, tbl$UNIPROT)
  syms <- map[acs]
  data.frame(
    input_ac = acs,
    symbol   = unname(syms),
    mapped   = !is.na(syms),
    stringsAsFactors = FALSE
  )
}

fbgn_to_symbol <- function(ids) {
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (!length(ids)) return(character(0))
  tbl <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys    = ids,
    keytype = "FLYBASE",
    columns = "SYMBOL"
  )
  tbl <- tbl[!duplicated(tbl$FLYBASE), ]
  map <- setNames(tbl$SYMBOL, tbl$FLYBASE)
  out <- map[ids]
  out[is.na(out)] <- ids[is.na(out)]
  out
}

# Build mitochondrial universe (symbols + UniProt ACs)
message("Building mitochondrial protein universe from GO annotations...")

mito_go_bp <- c("GO:0006119", "GO:0022900", "GO:0006099", "GO:0006635", "GO:0000422")
mito_go_cc <- c("GO:0005739", "GO:0005743", "GO:0005747", "GO:0005753", "GO:0005759")

.expand_offspring <- function(go_ids, offspr_env) {
  kids <- unique(unlist(mget(go_ids, envir = offspr_env, ifnotfound = NA)))
  unique(na.omit(c(go_ids, kids)))
}

all_mito_go <- unique(c(
  .expand_offspring(mito_go_bp, GOBPOFFSPRING),
  .expand_offspring(mito_go_cc, GOCCOFFSPRING)
))

go2fb <- AnnotationDbi::select(
  org.Dm.eg.db, keys = all_mito_go, keytype = "GO", columns = "FLYBASE"
)
mito_fbgn <- unique(na.omit(go2fb$FLYBASE))
mito_syms  <- unique(na.omit(unname(fbgn_to_symbol(mito_fbgn))))
message("  Symbol universe  : ", length(mito_syms), " gene symbols")
message("  Sample symbols   : ", paste(head(mito_syms, 10), collapse = ", "))

go2uni <- AnnotationDbi::select(
  org.Dm.eg.db, keys = all_mito_go, keytype = "GO", columns = "UNIPROT"
)
mito_acs <- unique(na.omit(go2uni$UNIPROT))
message("  UniProt AC universe: ", length(mito_acs), " ACs")

# Shared significance filter
.filter_sig <- function(df, id_col, label) {
  missing <- setdiff(c(id_col, PADJ_COL, LFC_COL), colnames(df))
  if (length(missing)) {
    stop(label, ": missing column(s): ", paste(missing, collapse = ", "),
         "\n  Available: ", paste(colnames(df), collapse = ", "))
  }
  df %>%
    filter(
      !is.na(.data[[PADJ_COL]]),
      .data[[PADJ_COL]] < PADJ_CUT,
      !is.na(.data[[LFC_COL]]),
      abs(.data[[LFC_COL]]) >= LFC_CUT
    ) %>%
    # Return the full filtered data frame (not just IDs) so logFC is available
    dplyr::select(all_of(c(id_col, PADJ_COL, LFC_COL))) %>%
    distinct(.data[[id_col]], .keep_all = TRUE)
}

# Diagnostic helper
.diag_trace <- function(raw_ids, parsed_acs = NULL, sym_df,
                         mito_syms, mito_acs = NULL, label, out_dir) {
  if (is.null(parsed_acs)) parsed_acs <- raw_ids

  in_mito_by_sym <- sym_df$symbol %in% mito_syms & !is.na(sym_df$symbol)
  in_mito_by_ac  <- if (!is.null(mito_acs)) parsed_acs %in% mito_acs else rep(FALSE, length(parsed_acs))
  in_mito        <- in_mito_by_sym | in_mito_by_ac

  trace <- data.frame(
    raw_id         = raw_ids,
    parsed_ac      = parsed_acs,
    symbol         = sym_df$symbol,
    mapped         = sym_df$mapped,
    in_mito_by_sym = in_mito_by_sym,
    in_mito_by_ac  = in_mito_by_ac,
    in_mito        = in_mito,
    stringsAsFactors = FALSE
  )

  n_total    <- nrow(trace)
  n_mapped   <- sum(trace$mapped, na.rm = TRUE)
  n_unmapped <- n_total - n_mapped
  n_mito     <- sum(trace$in_mito, na.rm = TRUE)
  n_sym_only <- sum(in_mito_by_sym & !in_mito_by_ac, na.rm = TRUE)
  n_ac_only  <- sum(!in_mito_by_sym & in_mito_by_ac, na.rm = TRUE)
  n_both_hit <- sum(in_mito_by_sym & in_mito_by_ac, na.rm = TRUE)

  message(label, " -- ID tracing summary:")
  message("    Total significant IDs    : ", n_total)
  message("    Mapped to a symbol       : ", n_mapped)
  message("    NOT mapped (no DB hit)   : ", n_unmapped)
  message("    In mito universe (total) : ", n_mito)
  message("      via symbol match only  : ", n_sym_only)
  message("      via AC match only      : ", n_ac_only)
  message("      via both               : ", n_both_hit)

  if (n_unmapped > 0)
    message("    Unmapped ACs: ",
            paste(trace$parsed_ac[!trace$mapped], collapse = ", "))
  if (n_mapped > 0 && n_mito == 0)
    message("    Mapped symbols (NOT mito): ",
            paste(na.omit(trace$symbol[trace$mapped]), collapse = ", "))
  if (n_mito > 0)
    message("    Mito-hit symbols: ",
            paste(na.omit(trace$symbol[trace$in_mito]), collapse = ", "))
  if (n_ac_only > 0)
    message("    Mito hits via AC only (no symbol): ",
            paste(trace$parsed_ac[!in_mito_by_sym & in_mito_by_ac], collapse = ", "))

  out_csv <- file.path(out_dir, "tables",
                       paste0("diag_ID_trace_", gsub("[^A-Za-z0-9]", "_", label), ".csv"))
  write.csv(trace, out_csv, row.names = FALSE)
  message("    Full trace saved -> ", out_csv)

  invisible(trace)
}

# Dataset loaders
# Each returns a data.frame: symbol | logFC (symbol falls back to UniProt AC if org.Dm.eg.db has no mapping)
load_dep_osd514 <- function(files, label) {
  dfs <- lapply(files, function(f) {
    if (!file.exists(f)) stop("File not found: ", f)
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  })
  df <- do.call(rbind, dfs)

  message(label, " -- columns: ", paste(colnames(df), collapse = ", "))
  message(label, " -- total rows: ", nrow(df))
  message(label, " -- sample protein_id: ", paste(head(df$protein_id, 5), collapse = ", "))
  message(label, " -- ", PADJ_COL, " range: [",
          round(min(df[[PADJ_COL]], na.rm=TRUE), 4), ", ",
          round(max(df[[PADJ_COL]], na.rm=TRUE), 4), "]")
  message(label, " -- ", LFC_COL, " range: [",
          round(min(df[[LFC_COL]], na.rm=TRUE), 3), ", ",
          round(max(df[[LFC_COL]], na.rm=TRUE), 3), "]")
  message(label, " -- rows passing padj<", PADJ_CUT, " only: ",
          sum(df[[PADJ_COL]] < PADJ_CUT, na.rm=TRUE))
  message(label, " -- rows passing |LFC|>=", LFC_CUT, " only: ",
          sum(abs(df[[LFC_COL]]) >= LFC_CUT, na.rm=TRUE))

  sig_df <- .filter_sig(df, "protein_id", label)
  message(label, " -- significant DEPs: ", nrow(sig_df))

  sym_df <- uniprot_to_symbol_df(sig_df$protein_id)
  trace  <- .diag_trace(sig_df$protein_id, NULL, sym_df, mito_syms, mito_acs, label, OUT_DIR)

  mito_rows <- which(trace$in_mito)
  out <- data.frame(
    symbol = ifelse(!is.na(trace$symbol[mito_rows]),
                    trace$symbol[mito_rows],
                    trace$parsed_ac[mito_rows]),
    logFC  = sig_df[[LFC_COL]][mito_rows],
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out$symbol), ]
  message(label, " -- mitochondrial DEPs: ", nrow(out))
  out
}

load_dep_tau <- function(files, label) {
  dfs <- lapply(files, function(f) {
    if (!file.exists(f)) stop("File not found: ", f)
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  })
  df <- do.call(rbind, dfs)

  message(label, " -- columns: ", paste(colnames(df), collapse = ", "))
  message(label, " -- total rows: ", nrow(df))
  message(label, " -- sample protein_id: ", paste(head(df$protein_id, 5), collapse = ", "))
  message(label, " -- ", PADJ_COL, " range: [",
          round(min(df[[PADJ_COL]], na.rm=TRUE), 4), ", ",
          round(max(df[[PADJ_COL]], na.rm=TRUE), 4), "]")
  message(label, " -- ", LFC_COL, " range: [",
          round(min(df[[LFC_COL]], na.rm=TRUE), 3), ", ",
          round(max(df[[LFC_COL]], na.rm=TRUE), 3), "]")
  message(label, " -- rows passing padj<", PADJ_CUT, " only: ",
          sum(df[[PADJ_COL]] < PADJ_CUT, na.rm=TRUE))
  message(label, " -- rows passing |LFC|>=", LFC_CUT, " only: ",
          sum(abs(df[[LFC_COL]]) >= LFC_CUT, na.rm=TRUE))

  sig_df <- .filter_sig(df, "protein_id", label)
  message(label, " -- significant DEPs: ", nrow(sig_df))

  acs    <- parse_uniprot_ac(sig_df$protein_id)
  sym_df <- uniprot_to_symbol_df(acs)
  trace  <- .diag_trace(sig_df$protein_id, acs, sym_df, mito_syms, mito_acs, label, OUT_DIR)

  mito_rows <- which(trace$in_mito)
  out <- data.frame(
    symbol = ifelse(!is.na(trace$symbol[mito_rows]),
                    trace$symbol[mito_rows],
                    trace$parsed_ac[mito_rows]),
    logFC  = sig_df[[LFC_COL]][mito_rows],
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out$symbol), ]
  message(label, " -- mitochondrial DEPs: ", nrow(out))
  out
}

# Extract mito DEP data frames
message("\n--- OSD-514 (Spaceflight) -- protein_id: bare UniProt ACs ---")
mito_osd514_df <- load_dep_osd514(OSD514_FILES, "OSD-514")

message("\n--- EmoryTau (Tau neurodegeneration) -- protein_id: FASTA headers ---")
mito_tau_df <- load_dep_tau(TAU_FILES, "EmoryTau Tau")

# Symbol vectors (for Venn)
mito_osd514 <- mito_osd514_df$symbol
mito_tau    <- mito_tau_df$symbol

# Save per-dataset protein lists (with logFC)
write.csv(mito_osd514_df,
          file.path(OUT_DIR, "tables", "mito_DEPs_OSD514.csv"), row.names = FALSE)
write.csv(mito_tau_df,
          file.path(OUT_DIR, "tables", "mito_DEPs_EmoryTau.csv"), row.names = FALSE)

# Build Venn list
venn_list <- list(
  "OSD-514\n(Spaceflight)" = mito_osd514,
  "EmoryTau\n(Tau)"        = mito_tau
)

# Draw Venn diagram
message("\nDrawing Venn diagram...")

venn         <- Venn(venn_list)
venn_regions <- process_region_data(venn)
message("venn_regions columns: ", paste(colnames(venn_regions), collapse = ", "))

p_venn <- ggVennDiagram(
  venn_list,
  label      = "count",
  label_size = 5,
  set_size   = 4.5
) +
  scale_fill_gradient(low = "#FEF3E2", high = "#D95F02") +
  scale_color_manual(values = c("#7B3F00", "#1B6CA8")) +

  labs(
    title    = "Mitochondrial DEPs \u2014 Cross-Dataset Overlap",
    subtitle = paste0(
      PADJ_COL, " < ", PADJ_CUT,
      "  |  |log2FC| \u2265 ", LFC_CUT,
      "  |  Mito protein universe: ", length(mito_syms), " genes"
    ),
    caption = "Numbers = mitochondrial DEP counts. See mito_DEP_venn_regions.csv for gene names."
  ) +

  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    plot.caption  = element_text(size = 8,  hjust = 0.5, color = "grey50"),
    legend.position = "right"
  )

out_png <- file.path(OUT_DIR, "figs", "mito_DEP_venn.png")
ggsave(out_png, p_venn, width = 10, height = 8, dpi = 300, bg = "white")
message("Venn diagram saved -> ", out_png)

# Companion table: gene names + logFC per Venn region
only_osd514 <- setdiff(mito_osd514, mito_tau)
only_tau    <- setdiff(mito_tau, mito_osd514)
both        <- intersect(mito_osd514, mito_tau)

max_len <- max(length(only_osd514), length(only_tau), length(both), 1)

region_table <- data.frame(
  OSD514_only   = c(only_osd514, rep(NA, max_len - length(only_osd514))),
  EmoryTau_only = c(only_tau,    rep(NA, max_len - length(only_tau))),
  Shared        = c(both,        rep(NA, max_len - length(both))),
  stringsAsFactors = FALSE
)

out_regions <- file.path(OUT_DIR, "tables", "mito_DEP_venn_regions.csv")
write.csv(region_table, out_regions, row.names = FALSE, na = "")
message("Region gene-name table saved -> ", out_regions)
message("  OSD-514 only  : ", length(only_osd514))
message("  EmoryTau only : ", length(only_tau))
message("  Shared        : ", length(both))

# Print overlapping proteins with direction in each dataset
message("\n=== Proteins shared between OSD-514 and EmoryTau (", length(both), " total) ===")

if (length(both) == 0) {
  message("  (none)")
} else {
  # Header
  message(sprintf("  %-20s  %-25s  %-25s", "Gene", "OSD-514", "EmoryTau (Tau)"))
  message(sprintf("  %-20s  %-25s  %-25s", "----", "-------", "--------------"))

  for (g in sort(both)) {
    lfc_osd <- mito_osd514_df$logFC[mito_osd514_df$symbol == g]
    lfc_tau <- mito_tau_df$logFC[mito_tau_df$symbol == g]

    dir_osd <- ifelse(lfc_osd > 0, "UP", "DOWN")
    dir_tau <- ifelse(lfc_tau > 0, "UP", "DOWN")

    message(sprintf("  %-20s  %-6s (logFC=%+.3f)    %-6s (logFC=%+.3f)",
                    g, dir_osd, lfc_osd, dir_tau, lfc_tau))
  }
}
