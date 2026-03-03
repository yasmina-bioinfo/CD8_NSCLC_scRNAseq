#!/usr/bin/env Rscript
# ============================================================
# 10_CD8_ModuleScore.R
# Module scores per CD8 state — PR vs SD
# Signatures derived from LUAD dataset for cross-dataset narrative
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

DefaultAssay(seu_cd8) <- "RNA"

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Palette
# -----------------------------
condition_colors <- c("PR" = "#4575B4", "SD" = "#D73027")

cluster_colors <- c(
  "CD8_Effector_GZMK"             = "#1f78b4",
  "CD8_Exhausted_Terminal"        = "#D73027",
  "CD8_Terminal_CX3CR1"           = "#FC8D59",
  "CD8_Proliferating"             = "#9970AB",
  "CD8_TRM_like"                  = "#8B0000",
  "CD8_Early_Activated_NR4A_high" = "#74c476",
  "CD8_Activated_HLAII_high"      = "#BCBD22",
  "CD8_IFN_Stress_Response"       = "#17BECF"
)

# -----------------------------
# Signatures — same as LUAD for cross-dataset narrative
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
clusters_keep <- c(
  "CD8_Effector_GZMK",
  "CD8_Exhausted_Terminal",
  "CD8_TRM_like",
  "CD8_Proliferating"
)

df <- seu_cd8@meta.data %>%
  select(cell_state, RECIST,
         Naive_score1, Effector_score1, TRM_score1) %>%
  filter(cell_state %in% clusters_keep) %>%
  mutate(
    cell_state = factor(cell_state, levels = clusters_keep),
    RECIST     = factor(RECIST, levels = c("PR", "SD"))
  )

# -----------------------------
# VlnPlot function
# -----------------------------
make_vln <- function(df, score_col, title, y_lab = "Module Score") {
  ggplot(df, aes(x = RECIST, y = .data[[score_col]],
                 fill = RECIST)) +
    geom_violin(trim = TRUE, alpha = 0.8, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3,
                 fill = "white", alpha = 0.7) +
    scale_fill_manual(values = condition_colors) +
    facet_wrap(~ cell_state, nrow = 1) +
    labs(title = title, y = y_lab, x = "") +
    theme_bw() +
    theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      strip.text      = element_text(size = 9,  face = "bold"),
      axis.text.x     = element_text(size = 10),
      axis.text.y     = element_text(size = 10),
      axis.title.y    = element_text(size = 11),
      legend.position = "none"
    )
}

# -----------------------------
# Three VlnPlots
# -----------------------------
p_naive    <- make_vln(df, "Naive_score1",    "Naive / Stem-like score")
p_effector <- make_vln(df, "Effector_score1", "Effector score")
p_TRM      <- make_vln(df, "TRM_score1",      "TRM Cytotoxic score")

# -----------------------------
# Combine
# -----------------------------
p_final <- p_naive / p_effector / p_TRM +
  plot_annotation(
    title = "CD8 functional state scores — PR vs SD",
    theme = theme(plot.title = element_text(size = 16, face = "bold",
                                            hjust = 0.5))
  )

# -----------------------------
# Save
# -----------------------------
ggsave(
  file.path(fig_dir, "10_CD8_ModuleScore.png"),
  p_final,
  width  = 16,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

cat("\nSaved -> Results/figures/10_CD8_ModuleScore.png\n")