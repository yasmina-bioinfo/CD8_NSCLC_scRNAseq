#!/usr/bin/env Rscript
# ============================================================
# 12b_CD8_CellChat_viz.R
# CellChat visualization — PR vs SD
# Anti-PD1 context
# ============================================================
library(CellChat)
library(ggplot2)
library(patchwork)

fig_dir <- "Results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load objects
# -----------------------------
cc_merged <- readRDS("Objects/12_CellChat_merged.rds")
cc_PR     <- readRDS("Objects/12_CellChat_PR.rds")
cc_SD     <- readRDS("Objects/12_CellChat_SD.rds")

# -----------------------------
# Palette
# -----------------------------
cluster_colors <- c(
  "CD8_Effector_GZMK"      = "#5BC8D4",
  "CD8_Exhausted_Terminal" = "#FD8D3C",
  "CD8_TRM_like"           = "#D73027",
  "CD8_Proliferating"      = "#9970AB",
  "T_Naive_Memory"         = "#1f78b4",
  "Cycling_Tcells"         = "#9970AB",
  "T_cells_Activated"      = "#FD8D3C",
  "B_cells"                = "#33a02c",
  "Plasma_cells"           = "#b2df8a",
  "TAM_like"               = "#e6550d",
  "TAM_like_MRC1"          = "#fd8d3c",
  "TAM_like_SPP1"          = "#fdae6b",
  "Monocytes_FCN1"         = "#fdd0a2",
  "Cycling_Myeloid"        = "#BCBD22",
  "Neutrophils"            = "#17BECF",
  "CAF_like"               = "#6a3d9a",
  "Endothelial_cells"      = "#b15928",
  "Mast_cells"             = "#fdbf6f",
  "Tumor_epithelial"       = "#a6cee3",
  "Tumor_epithelial_basal" = "#8b0000",
  "Tumor_epithelial_AT2"   = "#c51b8a",
  "Tumor_epithelial_EMT"   = "#de2d26",
  "Ciliated_epithelial"    = "#74c476"
)

# -----------------------------
# Plot 1 — Interaction count PR vs SD
# -----------------------------
p1 <- compareInteractions(
  cc_merged,
  show.legend = FALSE,
  group       = c(1, 2),
  measure     = "weight"
)

ggsave(file.path(fig_dir, "12a_CellChat_interactions_count.png"),
       p1, width = 6, height = 5, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "12a_CellChat_interactions_count.pdf"),
       p1, width = 6, height = 5, bg = "white")

# -----------------------------
# Plot 2 — Bubble plot PR vs SD
# -----------------------------
df_PR <- subsetCommunication(cc_PR,
                             sources.use = "CD8_Cytotoxic_Exhausted",
                             targets.use = c("TAM_like", "TAM_like_MRC1",
                                             "Tumor_epithelial", "Monocytes_FCN1"))
df_SD <- subsetCommunication(cc_SD,
                             sources.use = "CD8_Cytotoxic_Exhausted",
                             targets.use = c("TAM_like", "TAM_like_MRC1",
                                             "Tumor_epithelial", "Monocytes_FCN1"))

df_PR$condition <- "PR"
df_SD$condition <- "SD"
df_all <- rbind(df_PR, df_SD)
df_all <- df_all[df_all$prob > 0.05, ]

p_bubble <- ggplot(df_all,
                   aes(x = interaction_name_2,
                       y = target,
                       size = prob,
                       color = condition)) +
  geom_point(alpha = 1, stroke = 0.5) +
  scale_size_continuous(range = c(3, 12), name = "Prob") +
  scale_color_manual(values = c("PR" = "#4575B4", "SD" = "#D73027")) +
  facet_wrap(~ condition, ncol = 1) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 13, color = "#000000"),
    axis.text.y  = element_text(size = 13, color = "#000000"),
    strip.text   = element_text(size = 13, face = "bold", color = "#000000"),
    legend.title = element_text(size = 11, color = "#000000"),
    legend.text  = element_text(size = 10, color = "#000000"),
    plot.title   = element_text(size = 14, face = "bold", color = "#000000")
  ) +
  labs(x = "", y = "",
       title = "CD8_Exhausted_Terminal interactions — PR vs SD")

ggsave(file.path(fig_dir, "12c_CellChat_bubble_PRvsSD.png"),
       p_bubble, width = 14, height = 10, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "12c_CellChat_bubble_PRvsSD.pdf"),
       p_bubble, width = 14, height = 8, bg = "white")

cat("\nCellChat figures saved\n")