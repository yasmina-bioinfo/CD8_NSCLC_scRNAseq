#!/usr/bin/env Rscript

# ============================================================
# 11_Tcells_Markers_Top50.R
# Find markers per cluster for reclustered T cells (nLung/tLung)
# Outputs: full markers + top50 markers per cluster
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# -----------------------------
# Input / Output
# -----------------------------
in_rds <- "Objects/10_Tcells_nLung_tLung_reclustered.rds"

out_dir <- "Results/markers"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_all   <- file.path(out_dir, "11_Tcells_AllMarkers.csv")
out_top50 <- file.path(out_dir, "11_Tcells_Top50Markers.csv")

# -----------------------------
# Load object
# -----------------------------
seu_T <- readRDS(in_rds)

# Ensure identities are clusters
Idents(seu_T) <- "seurat_clusters"

# -----------------------------
# Find markers
# -----------------------------
markers_T <- FindAllMarkers(
  seu_T,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save full markers
write.table(markers_T, out_all, sep = "\t", quote = FALSE, row.names = FALSE)

# Top 50 per cluster
top50_T <- markers_T %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
  ungroup()

write.table(top50_T, out_top50, sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Minimal console validation
# -----------------------------
cat("\n=== 11_Tcells_Markers_Top50.R : VALIDATION ===\n")
cat("Object dims (features x cells): ", paste(dim(seu_T), collapse = " x "), "\n")
cat("Clusters:\n")
print(table(Idents(seu_T)))
cat("\nSaved:\n- ", out_all, "\n- ", out_top50, "\n")