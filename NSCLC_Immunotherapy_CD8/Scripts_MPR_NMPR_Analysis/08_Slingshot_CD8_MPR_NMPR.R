#!/usr/bin/env Rscript
# ============================================================
# 08_Slingshot_CD8_MPR_NMPR.R
# Pseudotime trajectory analysis on CD8 clusters
# GSE207422 — MPR vs NMPR (neoadjuvant anti-PD-1)
# Based on 18_Slingshot_CD8.R (LUAD dataset)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

# -----------------------------
# Load
# -----------------------------
seu_cd8 <- readRDS("Objects/07_CD8_annotated.rds")
DefaultAssay(seu_cd8) <- "RNA"

seu_cd8 <- subset(seu_cd8, subset = PathResponse %in% c("MPR", "NMPR"))

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Palette — 8 états CD8 GSE207422
# -----------------------------
cluster_colors <- c(
  "CD8_Activated_HLAII_high"       = "#A8D8E8",
  "CD8_Early_Activated_NR4A_high"  = "#FDBB84",
  "CD8_Effector_GZMK"              = "#FEE08B",
  "CD8_Exhausted_Terminal"         = "#D73027",
  "CD8_IFN_Stress_Response"        = "#7FCDBB",
  "CD8_Proliferating"              = "#CBC9E2",
  "CD8_Terminal_CX3CR1"            = "#4575B4",
  "CD8_TRM_like"                   = "#F46D43"
)

# -----------------------------
# Convertir en SCE
# -----------------------------
sce <- as.SingleCellExperiment(seu_cd8)

# -----------------------------
# Slingshot — racine = Early_Activated_NR4A_high
# Justification : état le moins différencié dans GSE207422
# (équivalent de Naive_CM dans LUAD)
# -----------------------------
sce <- slingshot(
  sce,
  clusterLabels = "cell_state",
  reducedDim    = "UMAP",
  start.clus    = "CD8_Early_Activated_NR4A_high"
)

# -----------------------------
# Extraire pseudotemps
# -----------------------------
pt <- slingPseudotime(sce)

if (ncol(pt) > 1) {
  seu_cd8$pseudotime <- rowMeans(pt, na.rm = TRUE)
} else {
  seu_cd8$pseudotime <- pt[, 1]
}

# -----------------------------
# Extraire les courbes
# -----------------------------
curves_list <- slingCurves(sce)
curves_df <- do.call(rbind, lapply(seq_along(curves_list), function(i) {
  df <- as.data.frame(curves_list[[i]]$s)
  colnames(df) <- c("umap_1", "umap_2")
  df$Lineage <- paste0("Lineage", i)
  df
}))

# -----------------------------
# UMAP dataframe
# -----------------------------
umap_df <- as.data.frame(Embeddings(seu_cd8, "umap"))
umap_df$cluster    <- seu_cd8$cell_state
umap_df$condition  <- seu_cd8$PathResponse  # MPR vs NMPR
umap_df$pseudotime <- seu_cd8$pseudotime

# -----------------------------
# Plot 1 — Trajectoire sur clusters annotés
# -----------------------------
p_traj <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 1.0, alpha = 0.7) +
  scale_color_manual(values = cluster_colors) +
  geom_path(data = curves_df,
            aes(x = umap_1, y = umap_2, group = Lineage),
            color = "black", linewidth = 1.0, alpha = 0.9,
            inherit.aes = FALSE) +
  labs(title = "CD8 Trajectory", color = "") +
  guides(color = guide_legend(
    override.aes = list(size = 6, alpha = 1, shape = 16)
  )) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text     = element_text(size = 10),
    legend.key.size = unit(8, "mm"),
    axis.title      = element_text(size = 11),
    axis.text       = element_text(size = 10)
  )

# -----------------------------
# Plot 2 — Pseudotime global
# -----------------------------
p_pseudo <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
  geom_point(size = 1.0, alpha = 0.7) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  geom_path(data = curves_df,
            aes(x = umap_1, y = umap_2, group = Lineage),
            color = "white", linewidth = 1.0, alpha = 0.9,
            inherit.aes = FALSE) +
  labs(title = "Pseudotime") +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 11),
    axis.title   = element_text(size = 11),
    axis.text    = element_text(size = 10)
  )

# -----------------------------
# Plot 3 — Pseudotime split MPR vs NMPR
# -----------------------------
umap_df$condition <- factor(umap_df$condition, levels = c("MPR", "NMPR"))

p_split <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
  geom_point(size = 1.0, alpha = 0.7) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  facet_wrap(~ condition) +
  geom_path(data = curves_df,
            aes(x = umap_1, y = umap_2, group = Lineage),
            color = "white", linewidth = 0.8, alpha = 0.9,
            inherit.aes = FALSE) +
  labs(title = "Pseudotime — MPR vs NMPR") +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text   = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 11),
    axis.title   = element_text(size = 11),
    axis.text    = element_text(size = 10)
  )

# -----------------------------
# Combine
# -----------------------------
p_final <- (p_traj | p_pseudo) / p_split +
  plot_annotation(
    title = "CD8 T cell pseudotime trajectory — MPR vs NMPR",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# -----------------------------
# Save
# -----------------------------
ggsave(
  file.path(fig_dir, "08_CD8_Slingshot_MPR_NMPR.png"),
  p_final,
  width  = 14,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

cat("\nSaved -> Results/figures/08_CD8_Slingshot_MPR_NMPR.png\n")
