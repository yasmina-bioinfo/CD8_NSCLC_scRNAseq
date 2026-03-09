#!/usr/bin/env Rscript

# ============================================================
# 07B_CD4_vs_CD8_top50_from_object.R
# - Load reclustered T cells object
# - Define CD4_like and CD8_like by cluster IDs
# - Exclude NK/B (clusters 14 and 15)
# - Compute top50 markers for each lineage (global CD4 vs CD8)
#
# Input : objects/05_Tcells_reclustered.rds
# Output:
#   - Results/markers/Tcells_CD4_vs_CD8_markers.tsv
#   - Results/markers/Tcells_top50_CD4_like.tsv
#   - Results/markers/Tcells_top50_CD8_like.tsv
#   - objects/06_Tcells_CD4CD8_only.rds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

# -----------------------------
# Paths
# -----------------------------
in_rds <- "objects/05_Tcells_reclustered.rds"
outdir <- "Results/markers"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_markers <- file.path(outdir, "Tcells_CD4_vs_CD8_markers.tsv")
out_top50_cd4 <- file.path(outdir, "Tcells_top50_CD4_like.tsv")
out_top50_cd8 <- file.path(outdir, "Tcells_top50_CD8_like.tsv")
out_rds <- "objects/06_Tcells_CD4CD8_only.rds"

# -----------------------------
# Load
# -----------------------------
seu_T <- readRDS(in_rds)
DefaultAssay(seu_T) <- "RNA"
Idents(seu_T) <- "seurat_clusters"

# -----------------------------
# Define cluster sets
# -----------------------------
cd4_clusters <- c("1","3","4","8","11")
cd8_clusters <- c("0","2","5","6","7","9","10","12","13")
exclude_clusters <- c("14","15")  # NK/γδ and B

# If clusters are "g0" style, auto-adapt
levs <- levels(Idents(seu_T))
if (any(grepl("^g\\d+$", levs))) {
  cd4_clusters <- paste0("g", cd4_clusters)
  cd8_clusters <- paste0("g", cd8_clusters)
  exclude_clusters <- paste0("g", exclude_clusters)
}

# -----------------------------
# Create lineage label
# -----------------------------
cl <- as.character(Idents(seu_T))
lineage <- rep("Other", length(cl))
lineage[cl %in% cd4_clusters] <- "CD4_like"
lineage[cl %in% cd8_clusters] <- "CD8_like"
lineage[cl %in% exclude_clusters] <- "Exclude"

seu_T$Tcell_lineage_simple <- factor(lineage, levels = c("CD4_like","CD8_like","Exclude","Other"))

# Keep only CD4_like + CD8_like (drop NK/B and any Other)
seu_L <- subset(seu_T, subset = Tcell_lineage_simple %in% c("CD4_like","CD8_like"))
Idents(seu_L) <- "Tcell_lineage_simple"

message("\nCounts:")
print(table(seu_L$Tcell_lineage_simple))

saveRDS(seu_L, out_rds)
message("Saved object -> ", out_rds)

# -----------------------------
# FindMarkers: CD4_like vs CD8_like
# -----------------------------
markers <- FindMarkers(
  seu_L,
  ident.1 = "CD4_like",
  ident.2 = "CD8_like",
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25
) %>%
  rownames_to_column("gene")

write_tsv(markers, out_markers)
message("Saved markers -> ", out_markers)

# Top50 CD4-like (positive logFC)
top50_cd4 <- markers %>%
  filter(avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50)

write_tsv(top50_cd4, out_top50_cd4)
message("Saved top50 CD4_like -> ", out_top50_cd4)

# Top50 CD8-like (negative logFC in this direction)
top50_cd8 <- markers %>%
  filter(avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%          # most negative first
  slice_head(n = 50)

write_tsv(top50_cd8, out_top50_cd8)
message("Saved top50 CD8_like -> ", out_top50_cd8)

message("\nDONE.")