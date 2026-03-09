#!/usr/bin/env Rscript
# ============================================================
# 07b_CellChat_MPR_NMPR_Visualization.R
# CellChat visualization — MPR vs NMPR
# Anti-PD1 neoadjuvant context
# ============================================================
library(CellChat)
library(ggplot2)
library(patchwork)

fig_dir <- "Results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load objects
# -----------------------------
cc_merged <- readRDS("Objects/07_CellChat_merged.rds")
cc_MPR    <- readRDS("Objects/07_CellChat_MPR.rds")
cc_NMPR   <- readRDS("Objects/07_CellChat_NMPR.rds")

# -----------------------------
# Palette
# -----------------------------
cluster_colors <- c(
  "CD8_Effector_GZMK"            = "#5BC8D4",
  "CD8_Exhausted_Terminal"       = "#FD8D3C",
  "CD8_TRM_like"                 = "#D73027",
  "CD8_Proliferating"            = "#9970AB",
  "CD8_IFN_Stress_Response"      = "#BCBD22",
  "CD8_Early_Activated_NR4A_high"= "#74c476",
  "CD8_Activated_HLAII_high"     = "#17BECF",
  "CD8_Terminal_CX3CR1"          = "#e6550d",
  "T_Naive_Memory"               = "#1f78b4",
  "Cycling_Tcells"               = "#9970AB",
  "T_cells_Activated"            = "#FD8D3C",
  "B_cells"                      = "#33a02c",
  "Plasma_cells"                 = "#b2df8a",
  "TAM_like"                     = "#e6550d",
  "TAM_like_MRC1"                = "#fd8d3c",
  "TAM_like_SPP1"                = "#fdae6b",
  "Monocytes_FCN1"               = "#fdd0a2",
  "Cycling_Myeloid"              = "#BCBD22",
  "Neutrophils"                  = "#17BECF",
  "CAF_like"                     = "#6a3d9a",
  "Endothelial_cells"            = "#b15928",
  "Mast_cells"                   = "#fdbf6f",
  "Tumor_epithelial"             = "#a6cee3",
  "Tumor_epithelial_basal"       = "#8b0000",
  "Tumor_epithelial_AT2"         = "#c51b8a",
  "Tumor_epithelial_EMT"         = "#de2d26",
  "Ciliated_epithelial"          = "#74c476"
)

# -----------------------------
# Plot 1 — Interaction count MPR vs NMPR
# -----------------------------
p1 <- compareInteractions(
  cc_merged,
  show.legend = FALSE,
  group       = c(1, 2),
  measure     = "weight"
)

ggsave(file.path(fig_dir, "07a_CellChat_interactions_count.png"),
       p1, width = 6, height = 5, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "07a_CellChat_interactions_count.pdf"),
       p1, width = 6, height = 5, bg = "white")

# -----------------------------
# Plot 2 — Bubble plot MPR vs NMPR
# CD8_Exhausted_Terminal as source
# -----------------------------
df_MPR <- subsetCommunication(cc_MPR,
                              sources.use = "CD8_Cytotoxic_Exhausted",
                              targets.use = c("TAM_like", "TAM_like_MRC1",
                                              "Tumor_epithelial", "Monocytes_FCN1"))
df_NMPR <- subsetCommunication(cc_NMPR,
                               sources.use = "CD8_Cytotoxic_Exhausted",
                               targets.use = c("TAM_like", "TAM_like_MRC1",
                                               "Tumor_epithelial", "Monocytes_FCN1"))

df_MPR$condition  <- "MPR"
df_NMPR$condition <- "NMPR"
df_all <- rbind(df_MPR, df_NMPR)
df_all <- df_all[df_all$prob > 0.05, ]

p_bubble <- ggplot(df_all,
                   aes(x = interaction_name_2,
                       y = target,
                       size = prob,
                       color = condition)) +
  geom_point(alpha = 1, stroke = 0.5) +
  scale_size_continuous(range = c(3, 12), name = "Prob") +
  scale_color_manual(values = c("MPR" = "#4575B4", "NMPR" = "#D73027")) +
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
       title = "CD8_Cytotoxic_Exhausted interactions — MPR vs NMPR")

ggsave(file.path(fig_dir, "07b_CellChat_bubble_MPRvsNMPR.png"),
       p_bubble, width = 14, height = 10, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "07b_CellChat_bubble_MPRvsNMPR.pdf"),
       p_bubble, width = 14, height = 8, bg = "white")

cat("\nCellChat MPR/NMPR figures saved\n")
