#!/usr/bin/env Rscript
# ============================================================
# 14_CD8_Slingshot.R
# Pseudotime trajectory analysis on CD8 clusters — PR vs SD
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

DefaultAssay(seu_cd8) <- "RNA"

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Palette — same 4 clusters as LUAD
# -----------------------------
cluster_colors <- c(
  "CD8_Effector_GZMK"      = "#92C5DE",
  "CD8_Exhausted_Terminal" = "#D73027",
  "CD8_TRM_like"           = "#D6A5A5",
  "CD8_Proliferating"      = "#CBC9E2"
)

# -----------------------------
# Subset to 4 clusters only
# -----------------------------
seu_cd8_sub <- subset(seu_cd8, 
                      cell_state %in% names(cluster_colors))
seu_cd8_sub$cell_state <- factor(seu_cd8_sub$cell_state, 
                                 levels = names(cluster_colors))

# -----------------------------
# Convert to SCE
# -----------------------------
sce <- as.SingleCellExperiment(seu_cd8_sub)

# -----------------------------
# Slingshot — root = CD8_Effector_GZMK (TCF7+, least differentiated)
# -----------------------------
sce <- slingshot(
  sce,
  clusterLabels = "cell_state",
  reducedDim    = "UMAP",
  start.clus    = "CD8_Effector_GZMK"
)

# -----------------------------
# Extract pseudotime
# -----------------------------
pt <- slingPseudotime(sce)

if (ncol(pt) > 1) {
  seu_cd8_sub$pseudotime <- rowMeans(pt, na.rm = TRUE)
} else {
  seu_cd8_sub$pseudotime <- pt[, 1]
}

# -----------------------------
# Extract curves
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
umap_df <- as.data.frame(Embeddings(seu_cd8_sub, "umap"))
umap_df$cluster    <- seu_cd8_sub$cell_state
umap_df$condition  <- seu_cd8_sub$RECIST
umap_df$pseudotime <- seu_cd8_sub$pseudotime

# -----------------------------
# Plot 1 — Trajectory on annotated clusters
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
    legend.text     = element_text(size = 11),
    legend.key.size = unit(8, "mm"),
    axis.title      = element_text(size = 11),
    axis.text       = element_text(size = 10)
  )

# -----------------------------
# Plot 2 — Global pseudotime
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
# Plot 3 — Pseudotime split PR vs SD
# -----------------------------
umap_df$condition <- factor(umap_df$condition, levels = c("PR", "SD"))

p_split <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
  geom_point(size = 1.0, alpha = 0.7) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  facet_wrap(~ condition) +
  geom_path(data = curves_df,
            aes(x = umap_1, y = umap_2, group = Lineage),
            color = "white", linewidth = 0.8, alpha = 0.9,
            inherit.aes = FALSE) +
  labs(title = "Pseudotime — PR vs SD") +
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
    title = "CD8 T cell pseudotime trajectory — PR vs SD",
    theme = theme(plot.title = element_text(size = 18, face = "bold", 
                                            hjust = 0.5))
  )

# -----------------------------
# Save
# -----------------------------
ggsave(
  file.path(fig_dir, "11_CD8_Slingshot.png"),
  p_final,
  width  = 14,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

cat("\nSaved -> Results/figures/11_CD8_Slingshot.png\n")