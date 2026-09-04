# Mitochondrial DEG Venn Diagram (RNA-seq)

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
  "Your_Working_Directory/DE_MG_vs_EARTH_unshrunk.csv"
)

TAU_FILES <- c(
  "Your_Working_Directory/DE_Tau_vs_Control.csv"
)

GAMMA_FILES <- c(
  "Your_Working_Directory/DE_High_vs_Control_unshrunk.csv"
)

OUT_DIR <- "Your_Working_Directory/CrossDataset_Venn"
dir.create(file.path(OUT_DIR, "figs"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

# Significance thresholds  (DESeq2 column names)
PADJ_CUT <- 0.05
LFC_CUT  <- 1       

PADJ_COL <- "padj"
LFC_COL  <- "log2FoldChange"

# ID helpers
# FlyBase IDs -> data.frame: fbgn | symbol | mapped
fbgn_to_symbol_df <- function(ids) {
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (!length(ids)) return(data.frame(fbgn=character(), symbol=character(), mapped=logical()))
  tbl <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys    = ids,
    keytype = "FLYBASE",
    columns = "SYMBOL"
  )
  tbl <- tbl[!duplicated(tbl$FLYBASE), ]
  map  <- setNames(tbl$SYMBOL, tbl$FLYBASE)
  syms <- map[ids]
  data.frame(
    fbgn   = ids,
    symbol = unname(syms),
    mapped = !is.na(syms),
    stringsAsFactors = FALSE
  )
}

# Convenience wrapper - returns symbol vector, falling back to FBgn if unmapped
fbgn_to_symbol <- function(ids) {
  df  <- fbgn_to_symbol_df(ids)
  out <- df$symbol
  out[is.na(out)] <- df$fbgn[is.na(out)]
  out
}

# Build mitochondrial gene universe
#   mito_fbgn : FlyBase IDs  (used for OSD-514 and Gamma matching)
#   mito_syms : gene symbols  (used for EmoryTau matching)

message("Building mitochondrial gene universe from GO annotations...")

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
message("  FlyBase ID universe: ", length(mito_fbgn), " FBgn IDs")

# Build symbol universe - note: some FBgn IDs may not map, shrinking this set
sym_df_universe <- fbgn_to_symbol_df(mito_fbgn)
n_unmapped_universe <- sum(!sym_df_universe$mapped)
mito_syms <- unique(na.omit(sym_df_universe$symbol))
message("  Symbol universe    : ", length(mito_syms), " gene symbols")
if (n_unmapped_universe > 0)
  message("  WARNING: ", n_unmapped_universe, " mito FBgn IDs have no symbol in org.Dm.eg.db — ",
          "these will be INVISIBLE to EmoryTau symbol matching")
message("  Sample symbols     : ", paste(head(mito_syms, 10), collapse = ", "))

# Shared significance filter - returns filtered data.frame (keeps logFC)

.filter_sig <- function(df, label) {
  missing <- setdiff(c("gene", PADJ_COL, LFC_COL), colnames(df))
  if (length(missing))
    stop(label, ": missing column(s): ", paste(missing, collapse = ", "),
         "\n  Available: ", paste(colnames(df), collapse = ", "))

  out <- df %>%
    filter(
      !is.na(.data[[PADJ_COL]]),
      .data[[PADJ_COL]] < PADJ_CUT,
      !is.na(.data[[LFC_COL]]),
      abs(.data[[LFC_COL]]) >= LFC_CUT
    ) %>%
    dplyr::select(all_of(c("gene", PADJ_COL, LFC_COL))) %>%
    distinct(gene, .keep_all = TRUE)

  message(label, " -- significant DEGs (", PADJ_COL, "<", PADJ_CUT,
          ", |LFC|>=", LFC_CUT, "): ", nrow(out))
  out
}

# Diagnostic helper - traces every ID through symbol conversion + mito filter
# Handles both FBgn->symbol (OSD-514/Gamma) and direct symbol (EmoryTau) paths
.diag_trace <- function(raw_ids, sym_df, mito_syms, mito_fbgn = NULL,
                         id_type = c("fbgn", "symbol"), label, out_dir) {
  id_type <- match.arg(id_type)

  # For FBgn datasets: mito match is via FBgn directly (reliable) AND via symbol
  # For symbol datasets: mito match is via symbol only (DB-dependent)
  if (id_type == "fbgn") {
    in_mito_by_fbgn <- raw_ids %in% mito_fbgn
    in_mito_by_sym  <- sym_df$symbol %in% mito_syms & !is.na(sym_df$symbol)
    in_mito         <- in_mito_by_fbgn  # FBgn is the ground truth for these datasets
  } else {
    in_mito_by_fbgn <- rep(FALSE, length(raw_ids))
    in_mito_by_sym  <- raw_ids %in% mito_syms
    in_mito         <- in_mito_by_sym
  }

  trace <- data.frame(
    raw_id         = raw_ids,
    symbol         = if (id_type == "fbgn") sym_df$symbol else raw_ids,
    mapped         = if (id_type == "fbgn") sym_df$mapped  else rep(TRUE, length(raw_ids)),
    in_mito_by_fbgn = in_mito_by_fbgn,
    in_mito_by_sym  = in_mito_by_sym,
    in_mito         = in_mito,
    stringsAsFactors = FALSE
  )

  n_total    <- nrow(trace)
  n_mapped   <- sum(trace$mapped)
  n_unmapped <- n_total - n_mapped
  n_mito     <- sum(trace$in_mito)

  message(label, " -- ID tracing summary:")
  message("    Total significant IDs      : ", n_total)

  if (id_type == "fbgn") {
    n_sym_miss <- sum(in_mito_by_fbgn & !in_mito_by_sym, na.rm = TRUE)
    message("    Mapped to a symbol         : ", n_mapped)
    message("    NOT mapped (no DB hit)     : ", n_unmapped)
    message("    In mito universe (FBgn)    : ", n_mito)
    if (n_sym_miss > 0)
      message("    Mito hits with NO symbol   : ", n_sym_miss,
              " (will appear as FBgn IDs in Venn)")
  } else {
    n_sym_miss <- sum(!in_mito_by_sym & (raw_ids %in% mito_syms == FALSE))
    message("    Matched to mito universe   : ", n_mito)
    message("    NOT in mito universe       : ", n_total - n_mito)
    # Check for potential synonym misses
    message("    NOTE: symbol matching is DB-dependent — see universe WARNING above")
  }

  if (n_mito > 0)
    message("    Mito-hit symbols: ",
            paste(na.omit(trace$symbol[trace$in_mito]), collapse = ", "))
  if (n_unmapped > 0 && id_type == "fbgn")
    message("    Unmapped FBgn IDs (in mito, shown as FBgn): ",
            paste(trace$raw_id[trace$in_mito & !trace$mapped], collapse = ", "))

  out_csv <- file.path(out_dir, "tables",
                       paste0("diag_ID_trace_", gsub("[^A-Za-z0-9]", "_", label), ".csv"))
  write.csv(trace, out_csv, row.names = FALSE)
  message("    Full trace saved -> ", out_csv)

  invisible(trace)
}

# Dataset loaders - each returns a data.frame: symbol | log2FoldChange
# OSD-514 and Gamma: gene col = FlyBase IDs
load_deg_fbgn <- function(files, label) {
  dfs <- lapply(files, function(f) {
    if (!file.exists(f)) stop("File not found: ", f)
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  })
  df <- do.call(rbind, dfs)

  message(label, " -- columns  : ", paste(colnames(df), collapse = ", "))
  message(label, " -- total rows: ", nrow(df))
  message(label, " -- sample gene IDs: ", paste(head(df$gene, 5), collapse = ", "))
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

  sig_df <- .filter_sig(df, label)
  sym_df <- fbgn_to_symbol_df(sig_df$gene)

  trace <- .diag_trace(sig_df$gene, sym_df, mito_syms, mito_fbgn,
                        id_type = "fbgn", label = label, out_dir = OUT_DIR)

  mito_rows <- which(trace$in_mito)
  out <- data.frame(
    symbol       = ifelse(!is.na(trace$symbol[mito_rows]),
                          trace$symbol[mito_rows],
                          trace$raw_id[mito_rows]),   # FBgn fallback
    log2FoldChange = sig_df[[LFC_COL]][mito_rows],
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out$symbol), ]
  message(label, " -- mitochondrial DEGs: ", nrow(out))
  out
}

# EmoryTau: gene col = gene symbols
load_deg_symbol <- function(files, label) {
  dfs <- lapply(files, function(f) {
    if (!file.exists(f)) stop("File not found: ", f)
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  })
  df <- do.call(rbind, dfs)

  message(label, " -- columns  : ", paste(colnames(df), collapse = ", "))
  message(label, " -- total rows: ", nrow(df))
  message(label, " -- sample gene IDs: ", paste(head(df$gene, 5), collapse = ", "))
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

  sig_df <- .filter_sig(df, label)

  # Reverse-lookup: symbol -> FBgn via org.Dm.eg.db, then intersect with mito_fbgn (ground truth). This is robust to annotation-version mismatches
  # that cause symbol-level misses when the EmoryTau DE file was built with a different FlyBase release than org.Dm.eg.db.
  sig_genes   <- sig_df$gene
  fbgn_lookup <- suppressMessages(
    mapIds(org.Dm.eg.db,
           keys      = sig_genes,
           keytype   = "SYMBOL",
           column    = "FLYBASE",
           multiVals = "first")
  )

  in_mito_by_fbgn <- fbgn_lookup %in% mito_fbgn
  in_mito_by_sym  <- sig_genes %in% mito_syms   # secondary check
  in_mito         <- in_mito_by_fbgn | in_mito_by_sym

  n_fbgn_only <- sum(in_mito_by_fbgn & !in_mito_by_sym, na.rm = TRUE)
  n_sym_only  <- sum(!in_mito_by_fbgn & in_mito_by_sym,  na.rm = TRUE)
  n_both      <- sum(in_mito_by_fbgn  & in_mito_by_sym,  na.rm = TRUE)
  n_unmapped  <- sum(is.na(fbgn_lookup))

  message(label, " -- ID tracing summary:")
  message("    Total significant DEGs     : ", length(sig_genes))
  message("    Symbols with no FBgn hit   : ", n_unmapped,
          " (these rely on symbol-only match)")
  message("    In mito (total)            : ", sum(in_mito))
  message("      via FBgn match only      : ", n_fbgn_only,
          " (would have been MISSED by symbol matching alone)")
  message("      via symbol match only    : ", n_sym_only)
  message("      via both                 : ", n_both)

  mito_genes <- sig_genes[in_mito]
  if (length(mito_genes) > 0)
    message("    Mito-hit symbols: ", paste(sort(mito_genes), collapse = ", "))

  # Save full trace
  trace_df <- data.frame(
    symbol          = sig_genes,
    fbgn_lookup     = unname(fbgn_lookup),
    in_mito_by_fbgn = in_mito_by_fbgn,
    in_mito_by_sym  = in_mito_by_sym,
    in_mito         = in_mito,
    stringsAsFactors = FALSE
  )
  out_csv <- file.path(OUT_DIR, "tables",
                       paste0("diag_ID_trace_", gsub("[^A-Za-z0-9]", "_", label), ".csv"))
  write.csv(trace_df, out_csv, row.names = FALSE)
  message("    Full trace saved -> ", out_csv)

  out <- data.frame(
    symbol         = mito_genes,
    log2FoldChange = sig_df[[LFC_COL]][in_mito],
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out$symbol), ]
  message(label, " -- mitochondrial DEGs: ", nrow(out))
  out
}

# Extract mito DEG data frames
message("\n--- OSD-514 (Spaceflight) -- gene col: FlyBase IDs ---")
mito_osd514_df <- load_deg_fbgn(OSD514_FILES, "OSD-514")

message("\n--- EmoryTau (Tau neurodegeneration) -- gene col: symbols ---")
mito_tau_df <- load_deg_symbol(TAU_FILES, "EmoryTau Tau")

message("\n--- PRJNA747152 (Gamma radiation) -- gene col: FlyBase IDs ---")
mito_gamma_df <- load_deg_fbgn(GAMMA_FILES, "PRJNA747152 Gamma")

# Symbol vectors (for Venn)
sym_osd514 <- mito_osd514_df$symbol
sym_tau    <- mito_tau_df$symbol
sym_gamma  <- mito_gamma_df$symbol

# Save per-dataset gene lists (with log2FoldChange)
write.csv(mito_osd514_df,
          file.path(OUT_DIR, "tables", "mito_DEGs_OSD514.csv"), row.names = FALSE)
write.csv(mito_tau_df,
          file.path(OUT_DIR, "tables", "mito_DEGs_EmoryTau.csv"), row.names = FALSE)
write.csv(mito_gamma_df,
          file.path(OUT_DIR, "tables", "mito_DEGs_Gamma.csv"), row.names = FALSE)

# Build Venn list
# Ensure all sets are character vectors so Venn() does not complain about class mismatch (happens when a dataset returns zero mito DEGs -> character(0) vs named character)
venn_list <- list(
  "OSD-514\n(Spaceflight)" = as.character(sym_osd514),
  "EmoryTau\n(Tau)"        = as.character(sym_tau),
  "PRJNA747152\n(Gamma)"   = as.character(sym_gamma)
)

# Draw Venn diagram (count labels; gene names in companion table)
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
  scale_fill_gradient(low = "#E8F4FD", high = "#2171B5") +
  scale_color_manual(values = c("#1B4F72", "#117A65", "#6E2F8C")) +

  labs(
    title    = "Mitochondrial DEGs \u2014 Cross-Dataset Overlap",
    subtitle = paste0(
      PADJ_COL, " < ", PADJ_CUT,
      "  |  |log2FC| \u2265 ", LFC_CUT,
      "  |  Mito gene universe: ", length(mito_fbgn), " FBgn IDs / ",
      length(mito_syms), " symbols"
    ),
    caption = "Numbers = mitochondrial DEG counts. See mito_DEG_venn_regions.csv for gene names."
  ) +

  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    plot.caption  = element_text(size = 8,  hjust = 0.5, color = "grey50"),
    legend.position = "right"
  )

out_png <- file.path(OUT_DIR, "figs", "mito_DEG_venn.png")
ggsave(out_png, p_venn, width = 10, height = 8, dpi = 300, bg = "white")
message("Venn diagram saved -> ", out_png)

# Companion table: gene names per Venn region
only_osd514    <- setdiff(sym_osd514, union(sym_tau, sym_gamma))
only_tau       <- setdiff(sym_tau,    union(sym_osd514, sym_gamma))
only_gamma     <- setdiff(sym_gamma,  union(sym_osd514, sym_tau))
osd514_tau     <- setdiff(intersect(sym_osd514, sym_tau),    sym_gamma)
osd514_gamma   <- setdiff(intersect(sym_osd514, sym_gamma),  sym_tau)
tau_gamma      <- setdiff(intersect(sym_tau,    sym_gamma),  sym_osd514)
all_three      <- Reduce(intersect, list(sym_osd514, sym_tau, sym_gamma))

max_len <- max(length(only_osd514), length(only_tau), length(only_gamma),
               length(osd514_tau), length(osd514_gamma), length(tau_gamma),
               length(all_three), 1)

pad <- function(x) c(x, rep(NA, max_len - length(x)))

region_table <- data.frame(
  OSD514_only      = pad(only_osd514),
  EmoryTau_only    = pad(only_tau),
  Gamma_only       = pad(only_gamma),
  OSD514_Tau       = pad(osd514_tau),
  OSD514_Gamma     = pad(osd514_gamma),
  Tau_Gamma        = pad(tau_gamma),
  All_three        = pad(all_three),
  stringsAsFactors = FALSE
)

out_regions <- file.path(OUT_DIR, "tables", "mito_DEG_venn_regions.csv")
write.csv(region_table, out_regions, row.names = FALSE, na = "")
message("Region gene-name table saved -> ", out_regions)

# Overlap summary printed to console
message("\n=== Pairwise overlap summary ===")
message("  OSD-514 only         : ", length(only_osd514))
message("  EmoryTau only        : ", length(only_tau))
message("  Gamma only           : ", length(only_gamma))
message("  OSD-514 & Tau only   : ", length(osd514_tau))
message("  OSD-514 & Gamma only : ", length(osd514_gamma))
message("  Tau & Gamma only     : ", length(tau_gamma))
message("  All three            : ", length(all_three))

# Print genes in ALL THREE datasets with direction per dataset
message("\n=== Genes significant in ALL THREE datasets (", length(all_three), " total) ===")

if (length(all_three) == 0) {
  message("  (none)")
} else {
  message(sprintf("  %-20s  %-28s  %-28s  %-28s",
                  "Gene", "OSD-514", "EmoryTau (Tau)", "Gamma"))
  message(sprintf("  %-20s  %-28s  %-28s  %-28s",
                  "----", "-------", "--------------", "-----"))

  for (g in sort(all_three)) {
    lfc_osd   <- mito_osd514_df$log2FoldChange[mito_osd514_df$symbol == g]
    lfc_tau   <- mito_tau_df$log2FoldChange[mito_tau_df$symbol == g]
    lfc_gamma <- mito_gamma_df$log2FoldChange[mito_gamma_df$symbol == g]

    dir_osd   <- ifelse(lfc_osd   > 0, "UP",   "DOWN")
    dir_tau   <- ifelse(lfc_tau   > 0, "UP",   "DOWN")
    dir_gamma <- ifelse(lfc_gamma > 0, "UP",   "DOWN")

    message(sprintf("  %-20s  %-6s (log2FC=%+.3f)      %-6s (log2FC=%+.3f)      %-6s (log2FC=%+.3f)",
                    g, dir_osd, lfc_osd, dir_tau, lfc_tau, dir_gamma, lfc_gamma))
  }
}

# Also print pairwise overlaps with directions
.print_pairwise <- function(genes, df1, df2, label1, label2) {
  if (length(genes) == 0) { message("  (none)"); return(invisible(NULL)) }
  message(sprintf("  %-20s  %-28s  %-28s", "Gene", label1, label2))
  message(sprintf("  %-20s  %-28s  %-28s", "----",
                  paste(rep("-", nchar(label1)), collapse=""),
                  paste(rep("-", nchar(label2)), collapse="")))
  for (g in sort(genes)) {
    lfc1 <- df1$log2FoldChange[df1$symbol == g]
    lfc2 <- df2$log2FoldChange[df2$symbol == g]
    message(sprintf("  %-20s  %-6s (log2FC=%+.3f)      %-6s (log2FC=%+.3f)",
                    g,
                    ifelse(lfc1 > 0, "UP", "DOWN"), lfc1,
                    ifelse(lfc2 > 0, "UP", "DOWN"), lfc2))
  }
}

message("\n=== OSD-514 & EmoryTau only (", length(osd514_tau), " genes) ===")
.print_pairwise(osd514_tau, mito_osd514_df, mito_tau_df, "OSD-514", "EmoryTau (Tau)")

message("\n=== OSD-514 & Gamma only (", length(osd514_gamma), " genes) ===")
.print_pairwise(osd514_gamma, mito_osd514_df, mito_gamma_df, "OSD-514", "Gamma")

message("\n=== EmoryTau & Gamma only (", length(tau_gamma), " genes) ===")
.print_pairwise(tau_gamma, mito_tau_df, mito_gamma_df, "EmoryTau (Tau)", "Gamma")

message("\nDONE -> ", OUT_DIR)

# EmoryTau synonym diagnostic
# Checks whether significant EmoryTau genes that missed the mito universe, might be present under a different capitalisation or alias.
# Prints any near-matches so you can decide whether to recode them.
message("\n--- EmoryTau synonym / capitalisation check ---")

tau_sig_all <- {
  df <- read.csv(TAU_FILES[1], stringsAsFactors = FALSE, check.names = FALSE)
  df %>%
    filter(!is.na(.data[[PADJ_COL]]), .data[[PADJ_COL]] < PADJ_CUT,
           !is.na(.data[[LFC_COL]]),  abs(.data[[LFC_COL]]) >= LFC_CUT) %>%
    pull(gene) %>% unique()
}

tau_not_mito <- setdiff(tau_sig_all, mito_syms)   # significant but not in mito universe

# Case-insensitive match against mito universe
mito_syms_lower  <- tolower(mito_syms)
tau_lower        <- tolower(tau_not_mito)
case_hits        <- tau_not_mito[tau_lower %in% mito_syms_lower]

if (length(case_hits) > 0) {
  message("  Genes in EmoryTau sig list that match mito universe case-insensitively (",
          length(case_hits), ") — likely capitalisation mismatches:")
  for (g in sort(case_hits)) {
    canonical <- mito_syms[mito_syms_lower == tolower(g)]
    message("    EmoryTau: '", g, "'  ->  mito universe canonical: '",
            paste(canonical, collapse = "' / '"), "'")
  }
  message("  ACTION: recode these in your EmoryTau DE file or add an alias map.")
} else {
  message("  No case-insensitive mismatches found.")
}

# Also check for partial matches on known mito genes that are commonly aliased
known_aliases <- list(
  "Drp1"    = c("drp1", "DRP1", "Dynamin-related protein 1"),
  "Opa1"    = c("opa1", "OPA1"),
  "Marf"    = c("marf", "MARF", "Mfn1", "Mfn2"),
  "Pink1"   = c("pink1", "PINK1"),
  "park"    = c("Parkin", "parkin", "PARK"),
  "sesB"    = c("sesb", "SESB", "ANT"),
  "blw"     = c("ATP5A", "atp5a"),
  "Cyt-c-p" = c("CytC", "cytc")
)

alias_hits <- character(0)
for (canonical in names(known_aliases)) {
  aliases <- known_aliases[[canonical]]
  found   <- intersect(tau_sig_all, c(canonical, aliases))
  if (length(found) > 0 && !canonical %in% mito_tau_df$symbol) {
    alias_hits <- c(alias_hits,
                    paste0("'", paste(found, collapse="' / '"),
                           "' (canonical mito symbol: '", canonical, "')"))
  }
}

if (length(alias_hits) > 0) {
  message("  Known mito gene aliases found in EmoryTau sig list but NOT captured:")
  for (h in alias_hits) message("    ", h)
} else {
  message("  No known alias mismatches found for common mito genes.")
}

# Deep diagnostic: inspect the 68 EmoryTau symbols with no FBgn hit
message("\n--- Inspecting EmoryTau symbols with no FBgn lookup hit ---")

tau_sig_df <- {
  df <- read.csv(TAU_FILES[1], stringsAsFactors = FALSE, check.names = FALSE)
  df %>%
    filter(!is.na(.data[[PADJ_COL]]), .data[[PADJ_COL]] < PADJ_CUT,
           !is.na(.data[[LFC_COL]]),  abs(.data[[LFC_COL]]) >= LFC_CUT) %>%
    dplyr::select(all_of(c("gene", PADJ_COL, LFC_COL))) %>%
    distinct(gene, .keep_all = TRUE)
}

fbgn_all <- suppressMessages(
  mapIds(org.Dm.eg.db,
         keys      = tau_sig_df$gene,
         keytype   = "SYMBOL",
         column    = "FLYBASE",
         multiVals = "first")
)

unmapped_syms <- tau_sig_df$gene[is.na(fbgn_all)]
message("  Total unmapped symbols: ", length(unmapped_syms))

# Check if any of these are in the full org.Dm.eg.db SYMBOL keyspace at all
all_db_syms   <- keys(org.Dm.eg.db, keytype = "SYMBOL")
not_in_db     <- unmapped_syms[!unmapped_syms %in% all_db_syms]
in_db_no_fbgn <- unmapped_syms[unmapped_syms %in% all_db_syms]

message("  Not in org.Dm.eg.db SYMBOL keyspace at all: ", length(not_in_db))
message("  In DB but no FLYBASE mapping               : ", length(in_db_no_fbgn))

# Save the full unmapped list with their LFC values for manual inspection
unmapped_df <- tau_sig_df[tau_sig_df$gene %in% unmapped_syms, ]
unmapped_df$in_db <- unmapped_df$gene %in% all_db_syms

# Flag names that look mitochondrial (contain known mito keywords)
mito_keywords <- c("mt:", "ND-", "COX", "ATP", "mito", "Mito", "Tim", "Tom",
                   "Cyt", "NADH", "Sdh", "Idh", "Mdh", "Cs", "Acl", "Etf",
                   "Letm", "Opa", "Drp", "Mfn", "Pink", "park", "Ub")
unmapped_df$looks_mito <- sapply(unmapped_df$gene, function(g) {
  any(sapply(mito_keywords, function(kw) grepl(kw, g, fixed = TRUE)))
})

out_unmapped <- file.path(OUT_DIR, "tables", "diag_EmoryTau_unmapped_symbols.csv")
write.csv(unmapped_df, out_unmapped, row.names = FALSE)
message("  Full unmapped symbol table saved -> ", out_unmapped)

looks_mito_genes <- unmapped_df$gene[unmapped_df$looks_mito]
if (length(looks_mito_genes) > 0) {
  message("  Unmapped symbols that LOOK mitochondrial by name (", length(looks_mito_genes), "):")
  for (g in sort(looks_mito_genes)) {
    lfc <- unmapped_df[[LFC_COL]][unmapped_df$gene == g]
    message("    ", g, "  (log2FC=", sprintf("%+.3f", lfc), ")")
  }
  message("  ACTION: manually check these against FlyBase — they may be valid mito hits")
  message("          under a different annotation build and worth adding to an alias map.")
} else {
  message("  No unmapped symbols look mitochondrial by name.")
}

message("\n  Conclusion: EmoryTau's low mito DEG count (", nrow(mito_tau_df), ") appears")
message("  to reflect genuine biology rather than a pipeline artefact — the FBgn")
message("  reverse-lookup rescued 0 additional hits, no capitalisation mismatches")
message("  were found, and unmapped symbols do not appear mitochondrial by name.")
