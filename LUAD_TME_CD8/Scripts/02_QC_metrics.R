#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scales)
})

# -----------------------------
# Paths
# -----------------------------
seurat_rds <- "objects/02_seurat_40k_with_QC.rds"
outdir <- "results/figures/QC"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load object
# -----------------------------
seu <- readRDS(seurat_rds)
DefaultAssay(seu) <- "RNA"

# Ensure percent.mt exists
if (!"percent.mt" %in% colnames(seu@meta.data)) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
}

md <- seu@meta.data

# ============================================================
# Theme (readable + infographic-friendly)
# ============================================================
theme_info <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text  = element_text(size = 12),
    plot.background = element_rect(fill = "white", color = NA)
  )

log10_labels <- trans_format("log10", math_format(10^.x))

# ============================================================
# Scatter: Genes vs UMI
# ============================================================
p_genes_umi <- ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.20, alpha = 0.12, color = "black") +
  scale_x_log10(labels = log10_labels) +
  labs(title = "Detected genes vs UMI", x = "UMI (log10)", y = "Detected genes") +
  theme_info

# ============================================================
# Scatter: Mito vs UMI
# ============================================================
p_mito_umi <- ggplot(md, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(size = 0.20, alpha = 0.12, color = "black") +
  scale_x_log10(labels = log10_labels) +
  labs(title = "Mitochondrial content vs UMI", x = "UMI (log10)", y = "Percent mitochondrial") +
  theme_info

# ============================================================
# Violin plots (single-group, no points)
# ============================================================
Idents(seu) <- factor(rep("Cells", ncol(seu)), levels = "Cells")

v_theme <- theme_info +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

v_genes <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0) +
  labs(title = "Detected genes per cell", y = "Detected genes") +
  v_theme

v_umi <- VlnPlot(seu, features = "nCount_RNA", pt.size = 0) +
  scale_y_log10() +
  labs(
    title = "UMI counts per cell",
    y = "UMI counts (log10)"
  ) +
  v_theme

v_mito <- VlnPlot(seu, features = "percent.mt", pt.size = 0) +
  labs(title = "Mitochondrial percentage", y = "Percent mitochondrial") +
  v_theme

# ============================================================
# Elbow plot
# ============================================================
p_elbow <- ElbowPlot(seu, ndims = 50) +
  labs(title = "PCA scree plot", x = "Principal component", y = "Standard deviation") +
  theme_info

# ============================================================
# Export helper (PDF + PNG, with cairo fallback)
# ============================================================
save_plot <- function(p, name, w = 6, h = 6, dpi = 450) {
  
  pdf_path <- file.path(outdir, paste0(name, ".pdf"))
  png_path <- file.path(outdir, paste0(name, ".png"))
  
  # ---- PDF (vector) ----
  # Try cairo first, fallback to base pdf
  ok_pdf <- FALSE
  try({
    if (capabilities("cairo")) {
      ggsave(pdf_path, plot = p, width = w, height = h, units = "in",
             device = cairo_pdf, bg = "white")
    } else {
      ggsave(pdf_path, plot = p, width = w, height = h, units = "in",
             device = "pdf", bg = "white")
    }
    ok_pdf <- file.exists(pdf_path)
  }, silent = TRUE)
  
  # ---- PNG (raster) ----
  ok_png <- FALSE
  try({
    ggsave(png_path, plot = p, width = w, height = h, units = "in",
           dpi = dpi, bg = "white")
    ok_png <- file.exists(png_path)
  }, silent = TRUE)
  
  if (!ok_pdf || !ok_png) {
    stop(
      "Export failed for: ", name,
      "\nPDF: ", pdf_path, " (exists=", ok_pdf, ")",
      "\nPNG: ", png_path, " (exists=", ok_png, ")"
    )
  }
  
  message("Saved: ", basename(pdf_path), " + ", basename(png_path))
}

# ============================================================
# Export figures
# ============================================================
save_plot(p_genes_umi, "QC_genes_vs_UMI")
save_plot(p_mito_umi,  "QC_mito_vs_UMI")
save_plot(v_genes,     "QC_violin_genes")
save_plot(v_umi,       "QC_violin_UMI")
save_plot(v_mito,      "QC_violin_mito")
save_plot(p_elbow,     "QC_elbow")

message("DONE: QC figures exported -> ", outdir)