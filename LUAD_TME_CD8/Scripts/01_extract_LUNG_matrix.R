#!/usr/bin/env Rscript

# ============================================================
# GSE131907 (LUAD) — Step 1
# Extract LUNG cells from raw UMI matrix and save as checkpoint
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
})

# ---- Paths ----
raw_path  <- "Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
out_path  <- "Data/GSE131907_Lung_only_UMI_matrix.rds"

# ---- Load raw matrix ----
raw_mat <- readRDS(raw_path)

# ---- Extract LUNG columns ----
lung_cells <- grep("LUNG", colnames(raw_mat), value = TRUE)
lung_mat <- raw_mat[, lung_cells]

# ---- Sanity checks ----
message("Raw matrix dim:  ", paste(dim(raw_mat), collapse = " x "))
message("LUNG matrix dim: ", paste(dim(lung_mat), collapse = " x "))
message("LUNG cells:      ", length(lung_cells))

# ---- Save checkpoint ----
saveRDS(lung_mat, out_path)
message("Saved -> ", out_path)