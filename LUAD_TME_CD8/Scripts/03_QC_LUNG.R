#!/usr/bin/env Rscript

# ============================================================
# GSE131907 — Step 3A
# QC metrics + QC plots (NO normalization, NO PCA)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

in_rds  <- "objects/02_seurat_LUNG_raw.rds"
out_rds <- "objects/03A_seurat_LUNG_QCmetrics.rds"

dir.create("Results", showWarnings = FALSE, recursive = TRUE)

# Plot style: readable when scaled down + white background on export
base_size <- 14
theme_set(theme_bw(base_size = base_size))
theme_update(
  plot.background  = element_rect(fill = "white", colour = NA),
  panel.background = element_rect(fill = "white", colour = NA)
)

seu <- readRDS(in_rds)

# QC metric
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Plots
p1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("QC: nCount_RNA vs percent.mt")
ggsave("Results/QC_scatter_nCount_vs_percentmt.png",
       plot = p1, width = 7, height = 6, units = "in", dpi = 450, bg = "white")

p2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("QC: nCount_RNA vs nFeature_RNA")
ggsave("Results/QC_scatter_nCount_vs_nFeature.png",
       plot = p2, width = 7, height = 6, units = "in", dpi = 450, bg = "white")

p3 <- FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.mt") +
  ggtitle("QC: nFeature_RNA vs percent.mt")
ggsave("Results/QC_scatter_nFeature_vs_percentmt.png",
       plot = p3, width = 7, height = 6, units = "in", dpi = 450, bg = "white")

saveRDS(seu, out_rds)
message("Saved -> ", out_rds)