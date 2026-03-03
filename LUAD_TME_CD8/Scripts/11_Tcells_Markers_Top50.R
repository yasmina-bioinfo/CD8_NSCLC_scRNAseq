#!/usr/bin/env Rscript

# ============================================================
# 11_TcellsClean_Markers_Top50.R
# Markers (Top 50) per cluster on CLEAN reclustered T cells
# Input: object produced by script 10 (clean UMAP/clusters)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# -----------------------------
# Input / Output
# -----------------------------
in_rds <- "Objects/10_Tcells_clean_reclustered.rds"   

out_dir <- "Results/markers"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_all   <- file.path(out_dir, "11_TcellsClean_AllMarkers.tsv")
out_top50 <- file.path(out_dir, "11_TcellsClean_Top50Markers.tsv")

# -----------------------------
# Load
# -----------------------------
seu_clean <- readRDS(in_rds)

# Make sure identities are clusters
Idents(seu_clean) <- "seurat_clusters"

# -----------------------------
# Find markers
# -----------------------------
markers <- FindAllMarkers(
  seu_clean,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.table(markers, out_all, sep = "\t", quote = FALSE, row.names = FALSE)

top50 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
  ungroup()

write.table(top50, out_top50, sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Validation (console only)
# -----------------------------
cat("\n=== 11_TcellsClean_Markers_Top50.R : VALIDATION ===\n")
cat("Object dims (features x cells): ", paste(dim(seu_clean), collapse = " x "), "\n")
cat("Clusters:\n")
print(table(Idents(seu_clean)))
cat("\nSaved:\n- ", out_all, "\n- ", out_top50, "\n")