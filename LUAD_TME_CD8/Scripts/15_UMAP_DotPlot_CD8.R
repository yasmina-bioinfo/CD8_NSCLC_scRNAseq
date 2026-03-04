#!/usr/bin/env Rscript

# ============================================================
# 15_CD8_UMAP_DotPlot_SplitCondition.R
# CD8-only: add explicit CD8_state labels, then export:
# - UMAP by CD8_state
# - UMAP split by Sample_Origin
# - DotPlot split by Sample_Origin (data-driven markers only)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# Input / Output
# -----------------------------
in_rds  <- "Objects/14_CD8_reclustered.rds"
out_rds <- "Objects/15_CD8_annotated_states.rds"

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load
# -----------------------------
seu_cd8 <- readRDS(in_rds)
DefaultAssay(seu_cd8) <- "RNA"

stopifnot("seurat_clusters" %in% colnames(seu_cd8@meta.data))
stopifnot("Sample_Origin" %in% colnames(seu_cd8@meta.data))

# Drop empty condition levels (optional but clean)
seu_cd8$Sample_Origin <- droplevels(seu_cd8$Sample_Origin)

# IMPORTANT: drop empty cluster levels if any
seu_cd8$seurat_clusters <- as.character(seu_cd8$seurat_clusters)

# Remove B-contaminated / doublet-like cluster
seu_cd8 <- subset(seu_cd8, subset = seurat_clusters != 3)
seu_cd8$seurat_clusters <- as.character(seu_cd8$seurat_clusters)

cat("\nClusters after removing B-like cluster 3:\n")
print(sort(unique(seu_cd8$seurat_clusters)))

# -----------------------------
# CD8 cluster -> state mapping (from your Top markers)
# -----------------------------
# Active clusters you reported after droplevels: 0,1,2,4


cat("\nClusters present (CD8):\n")
print(sort(unique(seu_cd8$seurat_clusters)))

# Mapping (edit here if you see an extra cluster)
cluster_map_cd8 <- c(
  "0" = "CD8_Naive_CM",
  "1" = "CD8_Effector_GZMK",
  "2" = "CD8_TRM_Cytotoxic",
  "4" = "CD8_Proliferating"
)

seu_cd8$CD8_state <- unname(cluster_map_cd8[seu_cd8$seurat_clusters])

# Fail fast with explicit info
if (any(is.na(seu_cd8$CD8_state))) {
  missing <- sort(unique(seu_cd8$seurat_clusters[is.na(seu_cd8$CD8_state)]))
  stop(paste("Unmapped seurat_clusters:", paste(missing, collapse = ", ")))
}

# -----------------------------
# Palette (distinct, readable)
# -----------------------------
cd8_cols <- c(
  "CD8_Naive_CM"       = "#A8D8E8",  
  "CD8_Effector_GZMK"  = "#FDBB84",  
  "CD8_TRM_Cytotoxic"  = "#D73027",  
  "CD8_Proliferating"  = "#CBC9E2"   
)

# -----------------------------
# Ensure UMAP exists (compute only if missing)
# -----------------------------
if (!"umap" %in% Reductions(seu_cd8)) {
  # Use same dims strategy as earlier (adjust only if you already decided otherwise)
  dims_use <- 1:15
  if (!"pca" %in% Reductions(seu_cd8)) {
    seu_cd8 <- NormalizeData(seu_cd8)
    seu_cd8 <- FindVariableFeatures(seu_cd8, nfeatures = 2000)
    seu_cd8 <- ScaleData(seu_cd8, features = VariableFeatures(seu_cd8))
    seu_cd8 <- RunPCA(seu_cd8, npcs = 30, features = VariableFeatures(seu_cd8))
  }
  seu_cd8 <- FindNeighbors(seu_cd8, dims = dims_use)
  # do NOT recluster here; keep existing clusters
  seu_cd8 <- RunUMAP(seu_cd8, dims = dims_use)
}

# -----------------------------
# UMAP by CD8_state
# -----------------------------
p_umap <- DimPlot(
  seu_cd8,
  reduction = "umap",
  group.by = "CD8_state",
  cols = cd8_cols,
  label = TRUE
) + theme_bw()

ggsave(file.path(fig_dir, "15_CD8_UMAP_states.png"),
       plot = p_umap, width = 7, height = 6, dpi = 450)
ggsave(file.path(fig_dir, "15_CD8_UMAP_states.pdf"),
       plot = p_umap, width = 7, height = 6)

# -----------------------------
# UMAP split by condition
# -----------------------------
p_umap_split <- DimPlot(
  seu_cd8,
  reduction = "umap",
  group.by = "CD8_state",
  split.by = "Sample_Origin",
  cols = cd8_cols,
  label = TRUE,
  ncol = 2
) + 
  theme_bw() +
  theme(
    legend.text = element_text(size = 13),
    legend.key.size = unit(6, "mm"),
    legend.spacing.y = unit(0.5, "cm"),
    strip.text = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

ggsave(file.path(fig_dir, "15_CD8_UMAP_states_splitByCondition.png"),
       plot = p_umap_split, width = 12, height = 6, dpi = 450)
ggsave(file.path(fig_dir, "15_CD8_UMAP_states_splitByCondition.pdf"),
       plot = p_umap_split, width = 12, height = 6)

# -----------------------------
# DotPlot split by condition (data-driven markers only)
# -----------------------------
markers_cd8_final <- c(
  "IL7R","TCF7",
  "GZMK","KLRG1",
  "GZMB","GNLY","ZNF683","CXCR6",
  "MKI67","UBE2C"
)

# Keep only genes that exist (strict dataset compliance)
markers_cd8_final <- markers_cd8_final[markers_cd8_final %in% rownames(seu_cd8)]
stopifnot(length(markers_cd8_final) > 0)

p_dot <- DotPlot(
  seu_cd8,
  features = markers_cd8_final,
  group.by = "CD8_state",
  split.by = "Sample_Origin"
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
  )

ggsave(file.path(fig_dir, "15_CD8_DotPlot_states_splitByCondition.png"),
       plot = p_dot, width = 12, height = 6, dpi = 450)
ggsave(file.path(fig_dir, "15_CD8_DotPlot_states_splitByCondition.pdf"),
       plot = p_dot, width = 12, height = 6)

# -----------------------------
# Validation + Save
# -----------------------------
cat("\n=== 15_CD8_UMAP_DotPlot_SplitCondition : VALIDATION ===\n")
cat("Cells per CD8_state:\n")
print(table(seu_cd8$CD8_state))

cat("\nCells per condition:\n")
print(table(seu_cd8$Sample_Origin))

cat("\nMarkers used in DotPlot:\n")
print(markers_cd8_final)

saveRDS(seu_cd8, out_rds)
cat("\nSaved -> ", out_rds, "\n", sep = "")