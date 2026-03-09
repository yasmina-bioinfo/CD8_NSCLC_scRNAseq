#!/usr/bin/env Rscript
# ============================================================
# 16_CD8_FeaturePlot.R
# Composite FeaturePlot split by condition
# Explicit CD8 remodeling narrative
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# Load annotated CD8 object
# -----------------------------
seu_cd8 <- readRDS("Objects/17_CD8_clean_annotated.rds")
DefaultAssay(seu_cd8) <- "RNA"

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Forcer l'ordre des niveaux de condition
# -----------------------------
seu_cd8$Sample_Origin <- factor(seu_cd8$Sample_Origin,
                                levels = c("nLung", "tLung"))

# -----------------------------
# Genes (validated in dataset)
# -----------------------------
genes_use <- c("TCF7", "GZMB", "MKI67", "TIGIT")
genes_use <- genes_use[genes_use %in% rownames(seu_cd8)]
stopifnot(length(genes_use) == 4)

# -----------------------------
# Subset par condition (une seule fois, hors lapply)
# -----------------------------
seu_n <- subset(seu_cd8, Sample_Origin == "nLung")
seu_t <- subset(seu_cd8, Sample_Origin == "tLung")

# -----------------------------
# Fonction helper pour un panel propre
# -----------------------------
make_panel <- function(seu_obj, gene, vmin, vmax, title_label, show_legend = TRUE) {
  p <- FeaturePlot(
    seu_obj,
    features   = gene,
    cols       = c("grey95", "#8B0000"),
    min.cutoff = vmin,
    max.cutoff = vmax,
    pt.size    = 1.2
  ) +
    ggtitle(title_label) +
    theme_bw() +
    theme(
      plot.title        = element_text(size = 18, face = "bold",
                                       hjust = 0.5, margin = margin(b = 2)),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 10),
      legend.title      = element_blank(),
      legend.text       = element_text(size = 10),
      legend.key.height = unit(12, "mm"),   # plus haute pour lisibilité
      legend.key.width  = unit(5, "mm"),
      plot.margin       = margin(4, 4, 4, 4)
    ) +
    # Forcer 3 breaks arrondis lisibles
    scale_color_gradient(
      low    = "grey95",
      high   = "#8B0000",
      limits = c(vmin, vmax),
      breaks = round(seq(vmin, vmax, length.out = 3), 2),  # <-- fix décimales
      labels = function(x) sprintf("%.2f", x)
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

# -----------------------------
# Generate split FeaturePlot
# -----------------------------
plot_list <- lapply(genes_use, function(g) {
  
  # Calcul min/max sur l'ensemble des cellules CD8
  # garantit que l'échelle est comparable entre nLung et tLung
  v    <- FetchData(seu_cd8, vars = g)[, 1]
  vmax <- unname(quantile(v,        0.99, na.rm = TRUE))
  vmin <- unname(quantile(v[v > 0], 0.05, na.rm = TRUE))
  
  p_n <- make_panel(seu_n, g, vmin, vmax, "nLung", show_legend = FALSE)
  p_t <- make_panel(seu_t, g, vmin, vmax, "tLung", show_legend = TRUE)
  
  # Combine nLung | tLung avec titre du gène au-dessus
  row_plot <- wrap_plots(p_n, p_t, ncol = 2)
  
  # Ajouter le titre du gène via un textGrob explicite
  title_grob <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = g, 
             size = 7,          # taille en unités ggplot (≈ 36pt)
             fontface = "bold",
             hjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  wrap_plots(title_grob, row_plot, 
             ncol = 1, 
             heights = c(0.06, 1))
})

# -----------------------------
# Combine vertically (4 rows)
# -----------------------------
p_final <- wrap_plots(plot_list, ncol = 1) &
  theme(plot.margin = margin(5, 5, 5, 5))

# -----------------------------
# Save
# -----------------------------
ggsave(
  file.path(fig_dir, "16_CD8_FeaturePlot_Composite_split.png"),
  p_final,
  width  = 12,
  height = 16,
  dpi    = 300,
  bg     = "white"
)

cat("\nSaved -> Results/figures/16_CD8_FeaturePlot_Composite_split.png\n")