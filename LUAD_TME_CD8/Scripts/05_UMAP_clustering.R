#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# -----------------------------
# Paths (INPUT = object already PCA-ready)
# -----------------------------
seurat_rds <- "objects/04_seurat_40k_norm_hvg_pca.rds"   
outdir <- "results/figures/UMAP"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load object
# -----------------------------
seu <- readRDS(seurat_rds)

# Safety checks
if (!"pca" %in% Reductions(seu)) {
  stop("ERROR: PCA not found in this Seurat object. Run PCA in Script 04 before Script 05.")
}

# -----------------------------
# Theme (same style as QC)
# -----------------------------
theme_info <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank()
  )

# -----------------------------
# Save function (PDF + PNG)
# -----------------------------
save_plot <- function(p, name, w = 6.5, h = 6.5) {
  # PDF vector
  ggsave(
    filename = file.path(outdir, paste0(name, ".pdf")),
    plot = p,
    width = w, height = h, units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  # PNG high-res
  ggsave(
    filename = file.path(outdir, paste0(name, ".png")),
    plot = p,
    width = w, height = h, units = "in",
    dpi = 450,
    bg = "white"
  )
}

# -----------------------------
# Parameters
# -----------------------------
dims_use <- 1:25
pt_size  <- 0.18   
pt_alpha <- 0.8   
do_raster <- TRUE  

# -----------------------------
# Neighbors + UMAP (if missing)
# -----------------------------
if (!"umap" %in% Reductions(seu)) {
  message("UMAP not found -> computing neighbors + UMAP...")
  seu <- FindNeighbors(seu, dims = dims_use, verbose = TRUE)
  seu <- RunUMAP(seu, dims = dims_use, verbose = TRUE)
} else {
 
  if (is.null(seu@graphs) || length(seu@graphs) == 0) {
    message("Neighbors graph not found -> computing neighbors...")
    seu <- FindNeighbors(seu, dims = dims_use, verbose = TRUE)
  }
}

# -----------------------------
# Clustering (two resolutions)
# -----------------------------
resolutions <- c(0.5)

for (res in resolutions) {
  res_key <- paste0("RNA_snn_res.", res)
  
  message("Clustering at resolution = ", res)
  seu <- FindClusters(seu, resolution = res, verbose = TRUE)
  
  # Make sure clusters are a factor (clean legend/labels)
  seu$seurat_clusters <- factor(seu@meta.data[[res_key]])
  
  ncl <- nlevels(seu$seurat_clusters)
  
  # Saturated, readable palette (better than glasbey here)
  cols <- hue_pal(l = 55, c = 100)(ncl)
}

  # -----------------------------
  # UMAP - no legend
  # -----------------------------
  
  ncl  <- length(levels(seu$seurat_clusters))
  levs <- levels(seu$seurat_clusters)
  
  cols <- setNames(viridis::turbo(ncl), levs)
  
  p_umap <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "seurat_clusters",
    pt.size = 0.25,
    shuffle = TRUE,
    raster = FALSE
  ) +
    scale_color_manual(values = cols) +
    ggtitle(paste0("UMAP (resolution = 0.5)")) +
    theme_info +
    NoLegend()
  
  p_umap_labeled <- LabelClusters(p_umap, id = "seurat_clusters", repel = TRUE, size = 4)
  save_plot(p_umap_labeled, "UMAP_clusters_res0.5_labeled")

# -----------------------------
# Save Seurat object with UMAP + clusterings
# -----------------------------
saveRDS(seu, file = "objects/05_seurat_with_umap_clusters.rds")

message("DONE: UMAP exported (PDF+PNG) -> ", outdir)
message("DONE: object saved -> objects/05_seurat_with_umap_clusters.rds")