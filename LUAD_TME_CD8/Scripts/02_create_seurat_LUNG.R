#!/usr/bin/env Rscript

# ============================================================
# GSE131907 — Step 2
# Create Seurat object from LUNG-only sparse UMI matrix
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

mat_path <- "Data/GSE131907_Lung_only_UMI_matrix_sparse.rds"
out_rds  <- "Data/02_seurat_LUNG_raw.rds"

lung_sparse <- readRDS(mat_path)

seu <- CreateSeuratObject(
  counts  = lung_sparse,
  project = "GSE131907_LUNG",
  assay   = "RNA"
)

rm(lung_sparse); gc()

message("Seurat dim: ", paste(dim(seu), collapse = " x "))
message("Cells:      ", ncol(seu))
message("Genes:      ", nrow(seu))

saveRDS(seu, out_rds)
message("Saved -> ", out_rds)