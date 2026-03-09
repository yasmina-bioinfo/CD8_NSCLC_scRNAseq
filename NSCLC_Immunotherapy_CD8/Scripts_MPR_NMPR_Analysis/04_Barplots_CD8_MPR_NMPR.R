# ============================================================
# 04_Barplot_CD8_MPR_NMPR.R
# CD8 cluster proportions by condition (MPR/NMPR)
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

seu_mpr <- readRDS("Objects/08_CD8_MPR_NMPR.rds")

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

# Proportions
df <- seu_mpr@meta.data %>%
  transmute(
    condition = PathResponse,
    cell_type = cell_state
  ) %>%
  filter(!is.na(condition), !is.na(cell_type))

df_prop <- df %>%
  count(condition, cell_type, name = "n") %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Plot
p_bar <- ggplot(df_prop, aes(x = condition, y = prop, fill = cell_type)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = cd8_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  ylab("CD8 cluster proportion") +
  theme_classic() +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 13, face = "bold"),
    axis.text.y      = element_text(size = 11),
    axis.title.y     = element_text(size = 12),
    legend.position  = "right",
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(0.5, "cm"),
    legend.spacing.y = unit(0.4, "cm"),
    legend.title     = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(fig_dir, "04_CD8_Barplot_proportions_MPR_NMPR.png"),
       p_bar, width = 10, height = 6, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "04_CD8_Barplot_proportions_MPR_NMPR.pdf"),
       p_bar, width = 10, height = 6, bg = "white")

cat("Done — CD8 barplot saved\n")