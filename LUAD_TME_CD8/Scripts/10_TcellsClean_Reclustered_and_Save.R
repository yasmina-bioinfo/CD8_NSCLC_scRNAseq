#!/usr/bin/env Rscript

# ============================================================
# 10_TcellsClean_Recluster_and_Save.R
# Recluster from scratch on CLEAN T-lineage cells and SAVE object
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

in_rds  <- "Objects/09_Tcells_clean_Tlineage.rds"
out_rds <- "Objects/10_Tcells_clean_reclustered.rds"

seu <- readRDS(in_rds)
DefaultAssay(seu) <- "RNA"

# Recompute from scratch
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30)

dims_use <- 1:15
seu <- FindNeighbors(seu, dims = dims_use)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = dims_use)

# Validation
cat("\n=== 10_TcellsClean_Recluster_and_Save : VALIDATION ===\n")
cat("Dims (features x cells): ", paste(dim(seu), collapse = " x "), "\n")
cat("Clusters:\n")
print(table(Idents(seu)))
cat("\nConditions:\n")
print(table(seu$Sample_Origin))

saveRDS(seu, out_rds)
cat("\nSaved -> ", out_rds, "\n")