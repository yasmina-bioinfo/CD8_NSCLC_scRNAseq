#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# -----------------------------
# Paths
# -----------------------------
seurat_rds <- "objects/05_seurat_with_umap_clusters.rds"
outdir <- "results/markers"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load object
# -----------------------------
seu <- readRDS(seurat_rds)

DefaultAssay(seu) <- "RNA"
Idents(seu) <- "seurat_clusters"

# -----------------------------
# Parameters
# -----------------------------
logfc_threshold <- 0.25
min_pct         <- 0.25
top_n           <- 50

# -----------------------------
# Find markers (all clusters)
# -----------------------------
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  logfc.threshold = logfc_threshold,
  min.pct = min_pct,
  test.use = "wilcox",
  verbose = TRUE
)

# -----------------------------
# Top 30 per cluster
# -----------------------------
top_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = top_n) %>%
  ungroup()

# -----------------------------
# Save full table
# -----------------------------
write.table(
  markers,
  file = file.path(outdir, "all_markers_clusters.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Save top markers
# -----------------------------
write.table(
  top_markers,
  file = file.path(outdir, "top50_markers_per_cluster.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("DONE: markers exported -> ", outdir)