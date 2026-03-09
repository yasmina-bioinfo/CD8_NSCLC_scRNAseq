#!/usr/bin/env Rscript
# ============================================================
# 09_CD8_FeaturePlot.R
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

DefaultAssay(seu_cd8) <- "RNA"

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Forcer l'ordre des niveaux
# -----------------------------
seu_cd8$RECIST <- factor(seu_cd8$RECIST, levels = c("PR", "SD"))

# -----------------------------
# Gènes — mêmes que LUAD
# -----------------------------
genes_use <- c("TCF7", "GZMB", "MKI67", "TIGIT")
genes_use <- genes_use[genes_use %in% rownames(seu_cd8)]
stopifnot(length(genes_use) == 4)

# -----------------------------
# Subset par condition
# -----------------------------
seu_PR <- subset(seu_cd8, RECIST == "PR")
seu_SD <- subset(seu_cd8, RECIST == "SD")

# -----------------------------
# Fonction helper
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
      legend.key.height = unit(8, "mm"),
      legend.key.width  = unit(5, "mm"),
      legend.margin     = margin(t = 0, b = 10)
    ) +
    scale_color_gradientn(
      colours = c("grey95", "#8B0000"),
      limits  = c(vmin, vmax),
      breaks  = c(vmin, vmax),
      labels  = c("0.00", sprintf("%.2f", vmax))
    )
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  return(p)
}

# -----------------------------
# Générer les panels
# -----------------------------
plot_list <- lapply(genes_use, function(g) {
  
  v    <- FetchData(seu_cd8, vars = g)[, 1]
  vmax <- unname(quantile(v, 0.99, na.rm = TRUE))
  vmin <- 0  # forcer à 0 comme LUAD
  
  p_PR <- make_panel(seu_PR, g, vmin, vmax, "PR", show_legend = FALSE)
  p_SD <- make_panel(seu_SD, g, vmin, vmax, "SD", show_legend = TRUE)
  
  # reste inchangé
  
  row_plot <- wrap_plots(p_PR, p_SD, ncol = 2)
  
  title_grob <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label    = g,
             size     = 7,
             fontface = "bold",
             hjust    = 0.5) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  wrap_plots(title_grob, row_plot,
             ncol    = 1,
             heights = c(0.06, 1))
})

# -----------------------------
# Assembler et sauvegarder
# -----------------------------
p_final <- wrap_plots(plot_list, ncol = 1) &
  theme(plot.margin = margin(5, 5, 5, 5))

ggsave(
  file.path(fig_dir, "CD8_FeaturePlot_split_PR_SD.png"),
  p_final,
  width  = 12,
  height = 16,
  dpi    = 300,
  bg     = "white"
)

cat("\nDone\n")