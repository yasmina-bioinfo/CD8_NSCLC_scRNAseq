#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

# -----------------------------
# Paths
# -----------------------------
input_rds  <- "objects/02_seurat_40k_with_QC.rds"
output_rds <- "objects/03_seurat_40k_filtered.rds"

# -----------------------------
# Load object
# -----------------------------
seu <- readRDS(input_rds)
DefaultAssay(seu) <- "RNA"

cat("Cells before filtering:", ncol(seu), "\n")

# ============================================================
# FILTERING CRITERIA
# nFeature: 500–5000
# percent.mt < 10%
# ============================================================

seu_filtered <- subset(
  seu,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 5000 &
    percent.mt < 10
)

cat("Cells after filtering:", ncol(seu_filtered), "\n")

# -----------------------------
# Save filtered object
# -----------------------------
saveRDS(seu_filtered, output_rds)

message("Filtering complete. Object saved -> ", output_rds)