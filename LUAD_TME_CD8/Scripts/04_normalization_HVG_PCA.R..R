#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scales)
})

# -----------------------------
# Paths
# -----------------------------
input_rds  <- "objects/03_seurat_40k_filtered.rds"
output_rds <- "objects/04_seurat_40k_norm_hvg_pca.rds"

figdir <- "results/figures/PCA"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load filtered object
# -----------------------------
seu <- readRDS(input_rds)
DefaultAssay(seu) <- "RNA"

cat("Cells (filtered):", ncol(seu), "\n")

# ============================================================
# Normalization -> HVG -> Scaling -> PCA
# ============================================================

seu <- NormalizeData(seu, verbose = FALSE)

seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

saveRDS(seu, "objects/04_seurat_40k_norm_hvg.rds")

seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)

seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)

# ============================================================
# Publication-quality theme
# ============================================================

theme_info <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    axis.line  = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# ============================================================
# Scree plot (post-filtering PCA)
# ============================================================

p_elbow <- ElbowPlot(seu, ndims = 50) +
  labs(
    title = "PCA scree plot (post-filtering)",
    x = "Principal component",
    y = "Standard deviation"
  ) +
  theme_info

# ============================================================
# Export function (PDF vector + PNG high-res)
# ============================================================

save_plot <- function(p, name, w = 6, h = 6, dpi = 450) {
  
  pdf_path <- file.path(figdir, paste0(name, ".pdf"))
  png_path <- file.path(figdir, paste0(name, ".png"))
  
  # PDF (vector)
  ggsave(
    filename = pdf_path,
    plot = p,
    width = w, height = h, units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  
  # PNG (high resolution)
  ggsave(
    filename = png_path,
    plot = p,
    width = w, height = h, units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  message("Saved: ", basename(pdf_path), " + ", basename(png_path))
}

# ============================================================
# Export figure
# ============================================================

save_plot(p_elbow, "PCA_scree_post_filtering")

# ============================================================
# PCA scatter (PC1 vs PC2)
# ============================================================

p_pca <- DimPlot(
  seu,
  reduction = "pca",
  dims = c(1, 2),
  pt.size = 0.4
) +
  labs(
    title = "PCA (PC1 vs PC2)",
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    axis.line  = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

save_plot(p_pca, "PCA_PC1_PC2")


# ============================================================
# Save processed object
# ============================================================

saveRDS(seu, output_rds)

message("DONE: normalized + PCA object saved -> ", output_rds)