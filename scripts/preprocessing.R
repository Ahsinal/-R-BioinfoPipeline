# =============================================================================
# 🧬 GENE ENRICHMENT PIPELINE - STEP 1: DATA PREPROCESSING
# =============================================================================

# PURPOSE:  →Extract differentially expressed genes (DEGs) from GEO datasets
#          Convert HGNC symbols → ENTREZ IDs for clusterProfiler analysis

# 1. load required packages
library(readr)      # Read TSV files
library(dplyr)      # Filter data  
library(clusterProfiler)  # Gene ID conversion (bitr function)
library(org.Hs.eg.db)     # Human genes annotation database

# 2. CREATE OUTPUT FOLDER
dir.create("data/processed/", showWarnings = FALSE)

# =============================================================================
# 🧪 CRITERIA FOR "IMPORTANT GENES":
# =============================================================================
# - `adj.P.Val` < 0.05  → Statistically significant (adjusted p-value)
# - `abs(logFC)` > 1     → Biologically meaningful (2-fold change up/down)

# 3. Process PARKINSON'S: Load → Filter → Convert
cat("📊Processing Parkinson's GSE1767...\n")
parkinson_raw <- read_tsv("data/raw/GSE1767.tsv")

# Display column names for verification
cat("📋 Available columns:", paste(names(parkinson_raw), collapse = ", "), "\n")

# Extract DEGs

# Filter: adj.P.val<0.05 (significant), |logFC|>1 (2x change)
parkinson_degs <- parkinson_raw %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  pull(Gene.symbol) %>% unique()

cat("✅ Parkinson's DEGs found:", length(parkinson_degs), "genes\n")

# Convert gene symbols to ENTREZ IDs (required for enrichment analysis)
parkinson_entrez <- bitr(parkinson_degs, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
cat("✅ Parkinson's ENTREZ conversion:", nrow(parkinson_entrez), "mapped\n\n")

# 4. Process HUNTINGTON'S: Load → Filter → Convert
cat("📊 Processing Huntington's GSE19587...\n")
huntington_raw <- read_tsv("data/raw/GSE19587.tsv")
cat("📋 Available columns:", paste(names(huntington_raw), collapse = ", "), "\n")

# Extract DEGs 
huntington_degs <- huntington_raw %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  pull(Gene.symbol) %>% unique()
cat("✅ Huntington's DEGs found:", length(huntington_degs), "genes\n")

# Convert to ENTREZ IDs
huntington_entrez <- bitr(huntington_degs, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)
cat("✅ Huntington's ENTREZ conversion:", nrow(huntington_entrez), "mapped\n\n")
# =============================================================================
# 5. SAVE PROCESSED RESULTS 
# =============================================================================
# Parkinson's results (symbols + ENTREZ IDs)
parkinson_results <- list(
  symbols = parkinson_degs,
  entrez_ids = parkinson_entrez$ENTREZID,
  total_degs = length(parkinson_degs),
  mapped_count = nrow(parkinson_entrez)
)

write_rds(parkinson_results, "data/processed/parkinson_degs.rds")


# Huntington's results
huntington_results <- list(
  symbols = huntington_degs,
  entrez_ids = huntington_entrez$ENTREZID,
  total_degs = length(huntington_degs),
  mapped_count = nrow(huntington_entrez)
)

write_rds(huntington_results, "data/processed/huntington_degs.rds")


cat("✅ Completed Preprocessing!\n")
cat("✅ Data SAVED! Ready for enrichment analysis\n")
cat("\n✅ Ready for GO/KEGG pathway analysis! 🎯\n")
