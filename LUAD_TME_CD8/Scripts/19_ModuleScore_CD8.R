#!/usr/bin/env Rscript
# ============================================================
# 18_CD8_ModuleScore.R
# Module scores per CD8 state — nLung vs tLung
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggpubr)
})

# -----------------------------
# Load
# -----------------------------
seu_cd8 <- readRDS("Objects/17_CD8_clean_annotated.rds")
DefaultAssay(seu_cd8) <- "RNA"

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Palette
# -----------------------------
condition_colors <- c("nLung" = "#4575B4", "tLung" = "#D73027")

cluster_colors <- c(
  "CD8_Naive_CM"       = "#5BC8D4",
  "CD8_Effector_GZMK"  = "#FD8D3C",
  "CD8_TRM_Cytotoxic"  = "#D73027",
  "CD8_Proliferating"  = "#9970AB"
)

# -----------------------------
# Signatures — filtrées pct.2 < 0.3
# -----------------------------
sig_naive <- list(c(
  "TCF7", "IL7R", "FOXP1", "LTB", "GPR183", "RORA", "CCR7"
))

sig_effector <- list(c(
  "GZMK", "GZMH", "GZMM", "NKG7", "CST7", "CCL4", "CCL3L3"
))

sig_TRM <- list(c(
  "GZMB", "GNLY", "ITGAE", "CXCR6", "ZNF683",
  "KLRC1", "TIGIT", "RGS1", "CTSW", "KLRD1"
))

# -----------------------------
# AddModuleScore
# -----------------------------
seu_cd8 <- AddModuleScore(seu_cd8, features = sig_naive,    name = "Naive_score")
seu_cd8 <- AddModuleScore(seu_cd8, features = sig_effector, name = "Effector_score")
seu_cd8 <- AddModuleScore(seu_cd8, features = sig_TRM,      name = "TRM_score")

# -----------------------------
# Metadata dataframe
# -----------------------------
df <- seu_cd8@meta.data %>%
  select(CD8_state, Sample_Origin,
         Naive_score1, Effector_score1, TRM_score1) %>%
  mutate(
    CD8_state      = factor(CD8_state, levels = names(cluster_colors)),
    Sample_Origin  = factor(Sample_Origin, levels = c("nLung", "tLung"))
  )

# -----------------------------
# Fonction VlnPlot propre
# -----------------------------
make_vln <- function(df, score_col, title, y_lab = "Module Score") {
  ggplot(df, aes(x = Sample_Origin, y = .data[[score_col]],
                 fill = Sample_Origin)) +
    geom_violin(trim = TRUE, alpha = 0.8, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3,
                 fill = "white", alpha = 0.7) +
    scale_fill_manual(values = condition_colors) +
    facet_wrap(~ CD8_state, nrow = 1) +
    labs(title = title, y = y_lab, x = "") +
    theme_bw() +
    theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      strip.text      = element_text(size = 11, face = "bold"),
      axis.text.x     = element_text(size = 10),
      axis.text.y     = element_text(size = 10),
      axis.title.y    = element_text(size = 11),
      legend.position = "none"
    )
}
# -----------------------------
# Trois VlnPlots
# -----------------------------
p_naive    <- make_vln(df, "Naive_score1",    "Naive / Stem-like score")
p_effector <- make_vln(df, "Effector_score1", "Effector score")
p_TRM      <- make_vln(df, "TRM_score1",      "TRM Cytotoxic score")

# -----------------------------
# Combine
# -----------------------------
p_final <- p_naive / p_effector / p_TRM +
  plot_annotation(
    title = "CD8 functional state scores — nLung vs tLung",
    theme = theme(plot.title = element_text(size = 16, face = "bold",
                                            hjust = 0.5))
  )

# -----------------------------
# Save
# -----------------------------
ggsave(
  file.path(fig_dir, "18_CD8_ModuleScore.png"),
  p_final,
  width  = 14,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

cat("\nSaved -> Results/figures/18_CD8_ModuleScore.png\n")