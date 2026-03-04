#!/usr/bin/env Rscript
# ============================================================
# 20_CellChat_Visualization_LUAD.R
# CellChat visualization — nLung vs tLung
# LUAD context
# ============================================================
library(CellChat)
library(ggplot2)
library(patchwork)

fig_dir <- "Results/figures"

# -----------------------------
# Load objects
# -----------------------------
cc_merged <- readRDS("Objects/19_CellChat_merged.rds")
cc_n      <- readRDS("Objects/19_CellChat_nLung.rds")
cc_t      <- readRDS("Objects/19_CellChat_tLung.rds")

# -----------------------------
# Palette
# -----------------------------
cluster_colors <- c(
  "CD8_Naive_CM"       = "#5BC8D4",
  "CD8_Effector_GZMK"  = "#FD8D3C",
  "CD8_TRM_Cytotoxic"  = "#D73027",
  "CD8_Proliferating"  = "#9970AB",
  "B lymphocytes"      = "#4DAF4A",
  "Epithelial cells"   = "#ADD8E6",
  "Myeloid cells"      = "#FF7F00",
  "Endothelial cells"  = "#8B4513",
  "Fibroblasts"        = "#7B2D8B",
  "MAST cells"         = "#FFFF99"
)

# -----------------------------
# Plot 1 — Interaction count nLung vs tLung
# -----------------------------
p1 <- compareInteractions(
  cc_merged,
  show.legend = FALSE,
  group       = c(1, 2),
  measure     = "weight"
)
ggsave(file.path(fig_dir, "19a_CellChat_interactions_count.png"),
       p1, width = 6, height = 5, dpi = 300, bg = "white")

# -----------------------------
# Plot 2 — Bubble plot nLung vs tLung
# CD8_TRM_Cytotoxic outgoing interactions
# -----------------------------
df_n <- subsetCommunication(cc_n,
                            sources.use = "CD8_TRM_Cytotoxic",
                            targets.use = c("Myeloid cells", "Epithelial cells",
                                            "B lymphocytes"))
df_t <- subsetCommunication(cc_t,
                            sources.use = "CD8_TRM_Cytotoxic",
                            targets.use = c("Myeloid cells", "Epithelial cells",
                                            "B lymphocytes"))

df_n$condition <- "nLung"
df_t$condition <- "tLung"

df_all <- rbind(df_n, df_t)

# Filter significant interactions

p_bubble <- ggplot(df_all,
                   aes(x = interaction_name_2,
                       y = target,
                       size = prob,
                       color = condition)) +
  geom_point(alpha = 1, stroke = 0.5) +
  scale_size_continuous(range = c(3, 12), name = "Prob") +
  scale_color_manual(values = c("nLung" = "#4575B4", "tLung" = "#D73027")) +
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
       title = "CD8_TRM_Cytotoxic interactions — nLung vs tLung")

ggsave(file.path(fig_dir, "19b_CellChat_bubble_nLung_tLung.png"),
       p_bubble, width = 14, height = 10, dpi = 450, bg = "white")

cat("\nCellChat LUAD figures saved\n")