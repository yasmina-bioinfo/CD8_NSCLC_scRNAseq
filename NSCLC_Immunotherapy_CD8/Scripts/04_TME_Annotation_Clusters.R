#!/usr/bin/env Rscript
# ============================================================
# Script 04 — Markers per cluster (TME annotation support)
# Output:
#  1) Results/markers/markers_all_clusters.tsv
#  2) Results/markers/top50_per_cluster.tsv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# -----------------------------
# Paths
# -----------------------------
obj_path <- "Objects/03_seu_umap_clustered.rds"  
out_dir  <- "results/markers"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# Load Seurat object
# -----------------------------
seu <- readRDS(obj_path)

# Assure-toi d'avoir des clusters
if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("seurat_clusters absent. Tu dois avoir fait FindNeighbors + FindClusters avant ce script.")
}
Idents(seu) <- "seurat_clusters"

DefaultAssay(seu) <- "RNA"

message("Cells: ", ncol(seu), " | Genes: ", nrow(seu), " | Clusters: ", length(levels(Idents(seu))))

# -----------------------------
# Find markers (RAM-friendly)
# -----------------------------
# NOTE: RAM friendly:
# - only.pos = TRUE
# - min.pct = 0.25
# - logfc.threshold = 0.25
# - test.use = "wilcox"
markers <- FindAllMarkers(
  object = seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Harmonize gene column,s name (Seurat can put rownames)
if (!"gene" %in% colnames(markers)) {
  markers <- markers %>% tibble::rownames_to_column("gene")
}

# Save
all_path <- file.path(out_dir, "markers_all_clusters.csv")
write_tsv(markers, all_path)
message("Saved: ", all_path)

# -----------------------------
# Top 50 markers per cluster
# -----------------------------
top_n <- 50
top50 <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = top_n) %>%
  ungroup()

top50_path <- file.path(out_dir, "top50_per_cluster.csv")
write_tsv(top50, top50_path)
message("Saved: ", top50_path)

message("DONE.")