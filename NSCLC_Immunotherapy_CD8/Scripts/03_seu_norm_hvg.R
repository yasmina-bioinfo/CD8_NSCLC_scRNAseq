#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Script 03 (low-RAM): Normalize -> HVG -> PCA -> UMAP -> Clustering
# Input : objects/02_seu_qc.rds
# Output: objects/03_seu_umap_clustered.rds + results/figures/*
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# -----------------------------
# Paths (EDIT THIS)
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy"
IN_OBJ   <- file.path(DATA_DIR, "objects", "02_seu_qc.rds")
OUT_OBJ  <- file.path(DATA_DIR, "objects")
OUT_FIG  <- file.path(DATA_DIR, "results", "figures")

dir.create(OUT_OBJ, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load QC object
# -----------------------------
seu <- readRDS(IN_OBJ)
DefaultAssay(seu) <- "RNA"

message("Loaded: ", IN_OBJ)
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

# -----------------------------
# 2) Normalize (low RAM)
# -----------------------------
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

# -----------------------------
# 3) HVG (keep standard 2000)
# -----------------------------
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Save checkpoint (optional but recommended)
saveRDS(seu, file.path(OUT_OBJ, "03a_seu_norm_hvg.rds"))
message("Saved: objects/03a_seu_norm_hvg.rds")

# -----------------------------
# 4) Scale only HVGs (IMPORTANT for RAM)
# -----------------------------
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)

# -----------------------------
# 5) PCA (30 dims enough)
# -----------------------------
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)

p_elbow <- ElbowPlot(seu, ndims = 30) + theme_bw()
ggsave(file.path(OUT_FIG, "03_ElbowPlot.png"), p_elbow, width = 5, height = 4, dpi = 300, bg = "white")

# -----------------------------
# 6) Neighbors + Clustering (start conservative)
# -----------------------------
dims_use <- 1:20   # adjust later after elbow
seu <- FindNeighbors(seu, dims = dims_use, verbose = FALSE)

# resolution: start 0.4; we can tune later
seu <- FindClusters(seu, resolution = 0.4, verbose = FALSE)

# -----------------------------
# 7) UMAP
# -----------------------------
seu <- RunUMAP(seu, dims = dims_use, verbose = FALSE)

# -----------------------------
# 8) Quick UMAP exports (by clusters + RECIST)
# -----------------------------
p_umap_clust <- DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE) + theme_bw()
ggsave(file.path(OUT_FIG, "03_UMAP_clusters.png"), p_umap_clust, width = 7, height = 5, dpi = 300, bg = "white")

p_umap_recist <- DimPlot(seu, reduction = "umap", group.by = "RECIST") + theme_bw()
ggsave(file.path(OUT_FIG, "03_UMAP_RECIST.png"), p_umap_recist, width = 7, height = 5, dpi = 300, bg = "white")

p_umap_split <- DimPlot(seu, reduction = "umap", split.by = "RECIST", ncol = 2) + theme_bw()
ggsave(file.path(OUT_FIG, "03_UMAP_split_RECIST.png"), p_umap_split, width = 12, height = 5, dpi = 300, bg = "white")

# -----------------------------
# 9) Save final object
# -----------------------------
saveRDS(seu, file.path(OUT_OBJ, "03_seu_umap_clustered.rds"))
message("Saved: objects/03_seu_umap_clustered.rds")
message("DONE Script 03")