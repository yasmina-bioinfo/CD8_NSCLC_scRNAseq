#!/usr/bin/env Rscript

# ============================================================
# 13_Tcells_Clean_Figures_UMAP_DotPlot.R
# Clean T-lineage object: recompute UMAP + export UMAP (global & split)
# and DotPlot validation
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# -----------------------------
# Input / Output
# -----------------------------
in_rds <- "Objects/12_Tcells_clean_Tlineage.rds"

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load
# -----------------------------
seu_clean <- readRDS(in_rds)
DefaultAssay(seu_clean) <- "RNA"
Idents(seu_clean) <- "seurat_clusters"

# -----------------------------
# Recompute UMAP (from scratch on clean object)
# -----------------------------
seu_clean <- NormalizeData(seu_clean)
seu_clean <- FindVariableFeatures(seu_clean, nfeatures = 2000)
seu_clean <- ScaleData(seu_clean)
seu_clean <- RunPCA(seu_clean, npcs = 30)

# Use PCs based on earlier elbow (~15). Keep consistent.
dims_use <- 1:15

seu_clean <- FindNeighbors(seu_clean, dims = dims_use)
seu_clean <- FindClusters(seu_clean, resolution = 0.5)
seu_clean <- RunUMAP(seu_clean, dims = dims_use)

# -----------------------------
# UMAP: clusters (global)
# -----------------------------
p_umap_clusters <- DimPlot(
  seu_clean,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) + theme_bw()

ggsave(file.path(fig_dir, "13_TcellsClean_UMAP_clusters.png"),
       plot = p_umap_clusters, width = 7, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "13_TcellsClean_UMAP_clusters.pdf"),
       plot = p_umap_clusters, width = 7, height = 6)

# -----------------------------
# UMAP: split by condition
# -----------------------------
p_umap_split <- DimPlot(
  seu_clean,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "Sample_Origin",
  label = TRUE,
  ncol = 2
) + theme_bw()

ggsave(file.path(fig_dir, "13_TcellsClean_UMAP_clusters_splitByCondition.png"),
       plot = p_umap_split, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "13_TcellsClean_UMAP_clusters_splitByCondition.pdf"),
       plot = p_umap_split, width = 12, height = 6)

# -----------------------------
# DotPlot: global validation markers
# -----------------------------
markers <- c(
  # General T
  "TRAC", "TRBC1", "TRBC2",
  # CD4 naive/memory
  "IL7R", "CCR7", "TCF7",
  # CD8 cytotoxic
  "CD8A", "CD8B", "GZMK", "NKG7", "PRF1", "GZMB",
  # Treg
  "FOXP3", "IL2RA", "CTLA4", "TIGIT"
)

p_dot <- DotPlot(seu_clean, features = markers) + RotatedAxis() + theme_bw()

ggsave(file.path(fig_dir, "13_TcellsClean_DotPlot_global.png"),
       plot = p_dot, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "13_TcellsClean_DotPlot_global.pdf"),
       plot = p_dot, width = 12, height = 6)

# -----------------------------
# Minimal console validation
# -----------------------------
cat("\n=== 13_Tcells_Clean_Figures_UMAP_DotPlot.R : VALIDATION ===\n")
cat("Dims (features x cells): ", paste(dim(seu_clean), collapse = " x "), "\n")
cat("Clusters:\n"); print(table(Idents(seu_clean)))
cat("\nConditions:\n"); print(table(seu_clean$Sample_Origin))
cat("\nFigures saved in: ", fig_dir, "\n")