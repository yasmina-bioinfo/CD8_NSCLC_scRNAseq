# =============================================================================
# Script 11 — Differentiation State Estimation with CytoTRACE
# Dataset  : Immunotherapy (GSE207422, Hu et al. 2023)
# Object   : 08_CD8_MPR_NMPR.rds (already loaded as seurat)
# Question : Are MPR CD8_Exhausted_Terminal cells less terminally differentiated
#            than NMPR — i.e. more plastic and potentially reactivable?
# Reference: Gulati G.S. et al. (2020). Science, 367(6476), 405-411.
#            doi:10.1126/science.aax0249
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------
# install.packages("https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz"

library(CytoTRACE)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)

# -----------------------------------------------------------------------------
# 1. LOAD OBJECT
# -----------------------------------------------------------------------------
seurat <- readRDS("C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/objects/08_CD8_MPR_NMPR.rds")

# -----------------------------------------------------------------------------
# 2. RUN CytoTRACE
# -----------------------------------------------------------------------------
mat_counts <- GetAssayData(seurat, assay = "RNA", layer = "counts")
mat_counts <- as.matrix(mat_counts)

results    <- CytoTRACE(mat_counts)
cyto_score <- results$CytoTRACE
seurat$CytoTRACE <- cyto_score[colnames(seurat)]

message("CytoTRACE scores computed for ", length(cyto_score), " cells")

# -----------------------------------------------------------------------------
# 3. SAVE OBJECT WITH SCORES
# -----------------------------------------------------------------------------
saveRDS(seurat, "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/Data/11_CD8_CytoTRACE.rds")
message("Saved: 11_CD8_CytoTRACE.rds")

# -----------------------------------------------------------------------------
# 4. PARAMETERS
# -----------------------------------------------------------------------------
fig_dir <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/Figures"

condition_colors <- c("MPR" = "#2166AC", "NMPR" = "#D6604D")

cd8_colors <- c(
  "CD8_Effector_GZMK"             = "#92C5DE",
  "CD8_Exhausted_Terminal"        = "#D73027",
  "CD8_Terminal_CX3CR1"           = "#FDBB84",
  "CD8_Proliferating"             = "#CBC9E2",
  "CD8_TRM_like"                  = "#D6A5A5",
  "CD8_Early_Activated_NR4A_high" = "#A1D99B",
  "CD8_Activated_HLAII_high"      = "#D9D97A",
  "CD8_IFN_Stress_Response"       = "#A8D8E8"
)

# -----------------------------------------------------------------------------
# 5. UMAP — CytoTRACE score with cluster labels
# -----------------------------------------------------------------------------

# Compute cluster centroids for labels
label_coords <- seurat@meta.data %>%
  mutate(
    UMAP1 = Embeddings(seurat, "umap")[, 1],
    UMAP2 = Embeddings(seurat, "umap")[, 2]
  ) %>%
  group_by(cell_state) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))

umap_cyto <- FeaturePlot(
  seurat,
  features   = "CytoTRACE",
  cols       = c("#2166AC", "lightgrey", "#D6604D"),
  order      = TRUE,
  min.cutoff = "q5",
  max.cutoff = "q95"
) +
  geom_label_repel(
    data          = label_coords,
    aes(x = UMAP1, y = UMAP2, label = cell_state),
    size          = 2.8,
    fontface      = "bold",
    color         = "black",
    fill          = alpha("white", 0.7),
    box.padding   = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size  = 0.3,
    inherit.aes   = FALSE
  ) +
  labs(
    title    = "CytoTRACE Differentiation Score — CD8 T cells",
    subtitle = "Blue = most differentiated | Red = most stem-like\n(Cluster annotation as in scRNA-seq landscape analysis)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle = element_text(size = 11, color = "#555555", face = "italic"),
    axis.text     = element_blank(),
    axis.ticks    = element_blank()
  )

ggsave(file.path(fig_dir, "11_UMAP_CytoTRACE_labeled.png"),
       umap_cyto, width = 8, height = 6, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "11_UMAP_CytoTRACE_labeled.pdf"),
       umap_cyto, width = 8, height = 6, bg = "white")

# -----------------------------------------------------------------------------
# 6. VIOLIN — CD8_Exhausted_Terminal MPR vs NMPR
# -----------------------------------------------------------------------------
exhausted <- seurat@meta.data %>%
  filter(cell_state == "CD8_Exhausted_Terminal", !is.na(CytoTRACE))

p_exhausted <- ggplot(exhausted, aes(x = PathResponse, y = CytoTRACE,
                                     fill = PathResponse)) +
  geom_violin(alpha = 0.8, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     size = 4, fontface = "italic", color = "#333333") +
  scale_fill_manual(values = condition_colors) +
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
    x        = "",
    y        = "CytoTRACE Score",
    title    = "Differentiation Potential — CD8_Exhausted_Terminal",
    subtitle = "MPR vs NMPR",
    caption  = "Wilcoxon test; cell-level comparison (patient-level validation pending)"
  )

ggsave(file.path(fig_dir, "11_Violin_CytoTRACE_Exhausted_MPR_NMPR.png"),
       p_exhausted, width = 6, height = 6, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "11_Violin_CytoTRACE_Exhausted_MPR_NMPR.pdf"),
       p_exhausted, width = 6, height = 6, bg = "white")

# -----------------------------------------------------------------------------
# 7. SESSION INFO
# -----------------------------------------------------------------------------
sink(file.path(fig_dir, "11_session_info.txt"))
sessionInfo()
sink()

message("\n=== Script 11 complete ===")
message("Outputs:")
message("  11_UMAP_CytoTRACE_labeled.png/pdf")
message("  11_Violin_CytoTRACE_Exhausted_MPR_NMPR.png/pdf")
message("  11_CD8_CytoTRACE.rds")