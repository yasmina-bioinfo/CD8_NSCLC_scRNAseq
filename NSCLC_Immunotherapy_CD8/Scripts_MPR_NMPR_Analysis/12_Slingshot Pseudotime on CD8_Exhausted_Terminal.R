# Script 12 — Slingshot Pseudotime on CD8_Exhausted_Terminal
# Subset: CD8_Exhausted_Terminal only
# Question: Do MPR and NMPR follow distinct differentiation trajectories
#           within the exhausted compartment?

library(Seurat)
library(slingshot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# -----------------------------------------------------------------------------
# 1. SUBSET CD8_Exhausted_Terminal
# -----------------------------------------------------------------------------

exhausted <- subset(seurat, cell_state == "CD8_Exhausted_Terminal")
message("Cells in subset: ", ncol(exhausted))
message("MPR: ", sum(exhausted$PathResponse == "MPR"))
message("NMPR: ", sum(exhausted$PathResponse == "NMPR"))

# -----------------------------------------------------------------------------
# 2. RE-CLUSTER WITHIN SUBSET
# -----------------------------------------------------------------------------
# Necessary to get meaningful sub-structure for trajectory

exhausted <- FindVariableFeatures(exhausted, nfeatures = 2000)
exhausted <- ScaleData(exhausted)
exhausted <- RunPCA(exhausted, npcs = 20)
exhausted <- RunUMAP(exhausted, dims = 1:15)
exhausted <- FindNeighbors(exhausted, dims = 1:15)
exhausted <- FindClusters(exhausted, resolution = 0.3)

# Quick check
DimPlot(exhausted, group.by = "seurat_clusters", label = TRUE)
DimPlot(exhausted, group.by = "PathResponse",
        cols = c("MPR" = "#2166AC", "NMPR" = "#D6604D"))

# -----------------------------------------------------------------------------
# 3. RUN SLINGSHOT
# -----------------------------------------------------------------------------

embed <- Embeddings(exhausted, "umap")
clusters <- exhausted$seurat_clusters

sds <- slingshot(
  data        = embed,
  clusterLabels = clusters,
  start.clus  = as.character(which.max(table(clusters)) - 1)
)

# Extract pseudotime
pt <- slingPseudotime(sds)
exhausted$pseudotime <- pt[, 1]

# -----------------------------------------------------------------------------
# 4. PLOT — Pseudotime colored by condition
# -----------------------------------------------------------------------------

df_plot <- data.frame(
  UMAP1      = embed[, 1],
  UMAP2      = embed[, 2],
  pseudotime = exhausted$pseudotime,
  condition  = exhausted$PathResponse
)

p_pseudo <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "C") +
  theme_bw(base_size = 13) +
  theme(
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    plot.title    = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle = element_text(size = 11, color = "#555555", face = "italic")
  ) +
  labs(
    title    = "Pseudotime — CD8_Exhausted_Terminal",
    subtitle = "Slingshot trajectory inference",
    color    = "Pseudotime"
  )

p_split <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "C") +
  facet_wrap(~ condition) +
  theme_bw(base_size = 13) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    strip.text       = element_text(size = 12, face = "bold", color = "#000000"),
    strip.background = element_rect(fill = "#F5F5F5"),
    plot.title       = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle    = element_text(size = 11, color = "#555555", face = "italic"),
    plot.caption     = element_text(size = 11, color = "#555555", face = "italic")
  ) +
  labs(
    title    = "Pseudotime — CD8_Exhausted_Terminal × Response",
    subtitle = "Slingshot trajectory inference",
    caption  = "Cell-level comparison; patient-level validation pending",
    color    = "Pseudotime"
  )

ggsave("12_Pseudotime_Exhausted_Terminal.png",
       p_pseudo, width = 8, height = 6, dpi = 450, bg = "white")
ggsave("12_Pseudotime_Exhausted_Terminal.pdf",
       p_pseudo, width = 8, height = 6, bg = "white")

ggsave("12_Pseudotime_Exhausted_Terminal_MPR_NMPR.png",
       p_split, width = 10, height = 6, dpi = 450, bg = "white")
ggsave("12_Pseudotime_Exhausted_Terminal_MPR_NMPR.pdf",
       p_split, width = 10, height = 6, bg = "white")

# -----------------------------------------------------------------------------
# 5. VIOLIN — Pseudotime distribution MPR vs NMPR
# -----------------------------------------------------------------------------

p_vln <- ggplot(df_plot, aes(x = condition, y = pseudotime,
                             fill = condition)) +
  geom_violin(alpha = 0.8, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("MPR" = "#2166AC", "NMPR" = "#D6604D")) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 12, color = "#000000"),
    axis.text.y      = element_text(size = 11, color = "#000000"),
    legend.position  = "none",
    plot.title       = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle    = element_text(size = 11, color = "#555555", face = "italic"),
    plot.caption     = element_text(size = 11, color = "#555555", face = "italic")
  ) +
  labs(
    title    = "Pseudotime Distribution — CD8_Exhausted_Terminal",
    subtitle = "MPR vs NMPR",
    caption  = "Wilcoxon test; cell-level comparison",
    x        = "",
    y        = "Pseudotime"
  ) +
  ggpubr::stat_compare_means(method = "wilcox.test",
                             label = "p.format", size = 4,
                             fontface = "italic", color = "#333333")

ggsave("12_Violin_Pseudotime_MPR_NMPR.png",
       p_vln, width = 6, height = 6, dpi = 450, bg = "white")
ggsave("12_Violin_Pseudotime_MPR_NMPR.pdf",
       p_vln, width = 6, height = 6, bg = "white")

message("\n=== Script 13 complete ===")
message("Outputs:")
message("  12_Pseudotime_Exhausted_Terminal.png/pdf")
message("  12_Pseudotime_Exhausted_Terminal_MPR_NMPR.png/pdf")
message("  12_Violin_Pseudotime_MPR_NMPR.png/pdf")