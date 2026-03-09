#!/usr/bin/env Rscript

# ============================================================
# 10_Reclustering_Tcells.R
# Reclustering from scratch on T cells (nLung + tLung)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# -----------------------------
# Input / Output
# -----------------------------
in_rds  <- "Objects/09_Tcells_nLung_tLung.rds"
out_rds <- "Objects/10_Tcells_nLung_tLung_reclustered.rds"

# -----------------------------
# Load
# -----------------------------
seu_T <- readRDS(in_rds)

# -----------------------------
# Recluster from scratch
# -----------------------------
DefaultAssay(seu_T) <- "RNA"

seu_T <- NormalizeData(seu_T)
seu_T <- FindVariableFeatures(seu_T, nfeatures = 2000)

seu_T <- ScaleData(seu_T)

seu_T <- RunPCA(seu_T, npcs = 30)

# NOTE: adjust dims after checking ElbowPlot interactively if needed
dims_use <- 1:15

seu_T <- FindNeighbors(seu_T, dims = dims_use)
seu_T <- FindClusters(seu_T, resolution = 0.5)

seu_T <- RunUMAP(seu_T, dims = dims_use)

# -----------------------------
# UMAP 
# -----------------------------

p_clusters <- DimPlot(
  seu_T,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) + theme_bw()

ggsave("Results/figures/10_Tcells_UMAP_clusters.png",
       plot = p_clusters, width = 7, height = 6, dpi = 300)

ggsave("Results/figures/10_Tcells_UMAP_clusters.pdf",
       plot = p_clusters, width = 7, height = 6)

# -----------------------------
# Minimal outputs for validation
# -----------------------------
cat("\n=== 10_Reclustering_Tcells.R : VALIDATION ===\n")
cat("Dims (features x cells): ", paste(dim(seu_T), collapse = " x "), "\n")
cat("Clusters:\n")
print(table(Idents(seu_T)))
cat("\nSample_Origin counts:\n")
print(table(seu_T$Sample_Origin))

# -----------------------------
# Save
# -----------------------------
dir.create("Objects", showWarnings = FALSE, recursive = TRUE)
saveRDS(seu_T, out_rds)
cat("\nSaved -> ", out_rds, "\n")