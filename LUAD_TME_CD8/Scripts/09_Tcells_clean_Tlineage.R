#!/usr/bin/env Rscript

# ============================================================
# 12_Clean_Tlineage.R
# Remove non-T clusters and evaluate T cell counts
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

# -----------------------------
# Input / Output
# -----------------------------
in_rds  <- "Objects/10_Tcells_nLung_tLung_reclustered.rds"
out_rds <- "Objects/10_Tcells_clean_reclustered.rds"

# -----------------------------
# Load object
# -----------------------------
seu_T <- readRDS(in_rds)

# Ensure identities are clusters
Idents(seu_T) <- "seurat_clusters"

# -----------------------------
# Keep only T-lineage clusters
# -----------------------------
T_clusters <- c("0", "1", "3", "5", "9")

seu_T_clean <- subset(seu_T, idents = T_clusters)

# -----------------------------
# Validation
# -----------------------------
cat("\n=== 12_Clean_Tlineage.R : VALIDATION ===\n")
cat("Dimensions (features x cells): ",
    paste(dim(seu_T_clean), collapse = " x "), "\n\n")

cat("Cells per cluster:\n")
print(table(Idents(seu_T_clean)))

cat("\nCells per condition:\n")
print(table(seu_T_clean$Sample_Origin))

# -----------------------------
# Save cleaned object
# -----------------------------
saveRDS(seu_T_clean, out_rds)

cat("\nSaved -> ", out_rds, "\n")