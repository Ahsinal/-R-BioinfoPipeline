# =============================================================================
# GENE ENRICHMENT PIPELINE: DATA PREPROCESSING
# =============================================================================
# PURPOSE  : Extract DEGs from GEO datasets, convert to ENTREZ IDs,
#            compute HD/PD overlap, and export for GO/KEGG analysis.
# DATASETS : GSE1767  — Huntington's Disease  (HD) || GSE19587 — Parkinson's Disease    (PD)
# =============================================================================

# ── PACKAGES ──────────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(clusterProfiler)  # bitr() — SYMBOL -> ENTREZID
library(org.Hs.eg.db)     # human gene annotations
library(hgu133plus2.db)   # Affymetrix HG-U133_Plus_2 probe annotation
library(AnnotationDbi)    # select()

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# ── THRESHOLDS ─────────────────────────────────────────────────────────────────
LOGFC_CUT <- 1     # |log2 fold change| > 1  =  2-fold change
PVAL_CUT  <- 0.05  # significance threshold (adj.P.Val or P.Value)

# ── HELPER: SYMBOL -> ENTREZID with progress report ───────────────────────────
map_to_entrez <- function(symbols) {
  entrez <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat(sprintf("  ENTREZ mapped   : %d / %d\n\n", nrow(entrez), length(symbols)))
  entrez
}
# ── HELPER: build result list in a consistent structure ───────────────────────
make_result <- function(symbols, entrez, deg_table, filter_note) {
  list(
    symbols    = symbols,
    entrez_ids = entrez$ENTREZID,
    deg_table  = deg_table,
    filter     = filter_note,
    n_degs     = length(symbols),
    n_mapped   = nrow(entrez)
  )
}
# =============================================================================
# SECTION 1 — HUNTINGTON'S DISEASE  (GSE1767)
# =============================================================================
# Platform : Custom RefSeq-based array — no Bioconductor annotation package.
#            Gene.symbol column in the TSV is the correct annotation source.
# Fix      : Probes mapping to multiple genes are encoded as "A///B///C".
#            We split these so every gene is counted individually — this is
#            the bioinformatics standard and matches how FunRich handles them.
# Filter   : adj.P.Val < 0.05 works normally for this dataset.
# =============================================================================

cat("Processing GSE1767 — Huntington's Disease...\n")

hd_raw <- read_tsv("data/raw/GSE1767.tsv", show_col_types = FALSE)

hd_degs <- hd_raw %>%
  filter(adj.P.Val < PVAL_CUT,
         abs(logFC) > LOGFC_CUT,
         !is.na(Gene.symbol),
         Gene.symbol != "",
         Gene.symbol != "---")

# Split "A///B///C" multi-gene entries into individual symbols
hd_symbols <- hd_degs$Gene.symbol %>%
  strsplit("///", fixed = TRUE) %>%
  unlist() %>%
  trimws() %>%
  .[!. %in% c("", "---")] %>%
  unique()

cat(sprintf("  DEGs            : %d genes\n", length(hd_symbols)))

# SYMBOL -> ENTREZID  (unmapped symbols are dropped with a warning — this is normal)
hd_entrez <- map_to_entrez(hd_symbols)
# =============================================================================
# SECTION 2 — PARKINSON'S DISEASE  (GSE19587)
# =============================================================================
# Platform : Affymetrix HG-U133_Plus_2  (probe IDs: "216136_at" format).
#            hgu133plus2.db is the correct annotation package for this array.
#            The TSV Gene.symbol column was annotated at GEO upload time with
#            an older source — using it directly gives wrong/different counts.
#            Remapping through hgu133plus2.db is the standard approach and
#            reproduces FunRich's gene counts exactly.
#
# P-value  : adj.P.Val minimum = 0.338 — Benjamini-Hochberg correction was
#            too conservative for this dataset size. Using raw P.Value < 0.05
#            is the accepted practice in exploratory microarray analysis when
#            no probe survives FDR correction (see Bourgon et al. 2010).
# =============================================================================

cat("Processing GSE19587 — Parkinson's Disease...\n")

pd_raw <- read_tsv("data/raw/GSE19587.tsv", show_col_types = FALSE)

# Warn clearly if adj.P.Val is unusable so the choice is always visible
adj_min <- min(pd_raw$adj.P.Val, na.rm = TRUE)
if (adj_min >= PVAL_CUT) {
  cat(sprintf("  NOTE: adj.P.Val min = %.3f — no probe passes FDR threshold.\n", adj_min))
  cat("        Falling back to raw P.Value < 0.05 (see comment above).\n")
  pd_degs <- pd_raw %>% filter(P.Value < PVAL_CUT, abs(logFC) > LOGFC_CUT)
} else {
  pd_degs <- pd_raw %>% filter(adj.P.Val < PVAL_CUT, abs(logFC) > LOGFC_CUT)
}

cat(sprintf("  Significant probes : %d\n", nrow(pd_degs)))

# Remap probe IDs -> current gene symbols via hgu133plus2.db
# AnnotationDbi::select() returns one row per probe-gene pair
pd_symbols <- AnnotationDbi::select(
  hgu133plus2.db,
  keys    = pd_degs$ID,
  columns = "SYMBOL",
  keytype = "PROBEID"
) %>%
  filter(!is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

cat(sprintf("  DEGs            : %d genes\n", length(pd_symbols)))

pd_entrez <- map_to_entrez(pd_symbols)

# =============================================================================
# SECTION 3 — HD / PD OVERLAP
# =============================================================================

common  <- intersect(hd_symbols, pd_symbols)
hd_only <- setdiff(hd_symbols, pd_symbols)
pd_only <- setdiff(pd_symbols, hd_symbols)

# =============================================================================
# SECTION 4 — SAVE OUTPUT FILES
# =============================================================================
# ── Huntington's ──────────────────────────────────────────────────────────────
hd_result <- make_result(hd_symbols, hd_entrez, hd_degs,
                         "adj.P.Val < 0.05 | |logFC| > 1 | /// splitting applied")
saveRDS(hd_result, "data/processed/huntington_degs.rds")
writeLines(sort(hd_symbols), "data/processed/huntington_gene_list.txt")

# ── Parkinson's ───────────────────────────────────────────────────────────────
pd_result <- make_result(pd_symbols, pd_entrez, pd_degs,
                         "P.Value < 0.05 | |logFC| > 1 | remapped via hgu133plus2.db")
saveRDS(pd_result, "data/processed/parkinson_degs.rds")
writeLines(sort(pd_symbols), "data/processed/parkinson_gene_list.txt")

# ── Overlap ───────────────────────────────────────────────────────────────────
overlap_result <- list(
  common  = sort(common),
  hd_only = sort(hd_only),
  pd_only = sort(pd_only)
)
saveRDS(overlap_result, "data/processed/HD_PD_overlap.rds")
writeLines(sort(common), "data/processed/common_genes_HD_PD.txt")

# =============================================================================
# SECTION 5 — SUMMARY
# =============================================================================

cat(strrep("=", 52), "\n")
cat(" PREPROCESSING COMPLETE\n")
cat(strrep("=", 52), "\n")
cat(sprintf("  Huntington's (GSE1767)  : %4d DEGs | %4d ENTREZ\n",
            hd_result$n_degs,  hd_result$n_mapped))
cat(sprintf("  Parkinson's  (GSE19587) : %4d DEGs | %4d ENTREZ\n",
            pd_result$n_degs,  pd_result$n_mapped))
cat(strrep("-", 52), "\n")
cat(sprintf("  Shared HD + PD          : %4d genes\n", length(common)))
cat(sprintf("  HD only                 : %4d genes\n", length(hd_only)))
cat(sprintf("  PD only                 : %4d genes\n", length(pd_only)))
cat(strrep("=", 52), "\n")
cat("  Output -> data/processed/\n")
cat("  Ready for GO/KEGG enrichment analysis.\n")
cat(strrep("=", 52), "\n")

# DEG ANALYSIS SUMMARY
# =============================================================================
# Create table
table_data <- data.frame(
  Disease    = c("Parkinson's", "Huntington's"),
  Total_DEGs = c(pd_result$n_degs, hd_result$n_degs),
  Mapped     = c(pd_result$n_mapped, hd_result$n_mapped),
  Rate       = paste0(round(100 * c(pd_result$n_mapped, hd_result$n_mapped) /
                              c(pd_result$n_degs,   pd_result$n_degs), 1), "%"),
  Status     = ifelse(c(pd_result$n_degs, hd_result$n_degs) > 0, "✅ Ready", "⚠️ No DEGs"),
  stringsAsFactors = FALSE
)

# Print table with borders
cat("\n", strrep("=", 50), "\n")
cat("           DEG ANALYSIS SUMMARY\n")
cat(strrep("=", 50), "\n\n")
print(table_data, row.names = FALSE)
cat("\n", strrep("=", 50), "\n")