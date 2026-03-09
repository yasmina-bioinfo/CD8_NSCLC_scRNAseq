#!/usr/bin/env Rscript
# ============================================================
# 02_UMAP_MPR_NMPR.R
# UMAP colored by cluster and by condition (MPR/NMPR)
# ============================================================
library(Seurat)
library(ggplot2)

fig_dir <- "Results/figures/MPR_NMPR"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

seu_mpr <- readRDS("Objects/08_CD8_MPR_NMPR.rds")

# Colors cluster
cluster_colors <- c(
  "CD8_Effector_GZMK"             = "#92C5DE",
  "CD8_Exhausted_Terminal"        = "#D73027",
  "CD8_Terminal_CX3CR1"           = "#FDBB84",
  "CD8_Proliferating"             = "#CBC9E2",
  "CD8_TRM_like"                  = "#D6A5A5",
  "CD8_Early_Activated_NR4A_high" = "#A1D99B",
  "CD8_Activated_HLAII_high"      = "#D9D97A",
  "CD8_IFN_Stress_Response"       = "#A8D8E8"
)

# UMAP by cluster
p1 <- DimPlot(seu_mpr, group.by = "cell_state",
              cols = cluster_colors, pt.size = 0.3) +
  ggtitle("CD8 clusters — MPR/NMPR") +
  theme(plot.title = element_text(size = 13, face = "bold"))

# UMAP by condition
p2 <- DimPlot(seu_mpr, group.by = "PathResponse",
              cols = c("MPR" = "#4575B4", "NMPR" = "#D73027"),
              pt.size = 0.3) +
  ggtitle("Condition — MPR vs NMPR") +
  theme(plot.title = element_text(size = 13, face = "bold"))

# UMAP split by condition with cluster annotations
p3 <- DimPlot(seu_mpr, group.by = "cell_state",
              split.by = "PathResponse",
              cols = cluster_colors, pt.size = 0.3,
              label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle("CD8 clusters — split by condition") +
  theme(plot.title = element_text(size = 13, face = "bold"),
        strip.text  = element_text(size = 12, face = "bold"))

# Save
ggsave(file.path(fig_dir, "02_UMAP_clusters_MPR_NMPR.png"),
       p1, width = 8, height = 6, dpi = 450, bg = "white")

ggsave(file.path(fig_dir, "02_UMAP_condition_MPR_NMPR.png"),
       p2, width = 7, height = 6, dpi = 450, bg = "white")

ggsave(file.path(fig_dir, "02_UMAP_split_annotated_MPR_NMPR.png"),
       p3, width = 14, height = 6, dpi = 450, bg = "white")

cat("Done — UMAPs saved\n")