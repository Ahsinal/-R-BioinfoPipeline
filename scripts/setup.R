# Initial project setup for Parkinson's/Huntington's enrichment analysis
cat("Installing required packages...")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "gprofiler2", "readr", "dplyr", "ggplot2"))
renv::snapshot()  # Lock versions

cat("Ready!\n")
