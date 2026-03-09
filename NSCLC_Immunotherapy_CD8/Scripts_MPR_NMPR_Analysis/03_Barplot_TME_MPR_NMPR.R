# ============================================================
# 03_Barplot_TME_MPR_NMPR.R
# Cell type proportions by condition (MPR/NMPR)
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

seu_tme_mpr <- readRDS("Objects/04_TME_MPR_NMPR.rds")

immuno_colors <- c(
  "T_Naive_Memory"           = "#1f78b4",
  "CD8_Cytotoxic_Exhausted"  = "#D73027",
  "Cycling_Tcells"           = "#9970AB",
  "T_cells_Activated"        = "#FD8D3C",
  "TAM_like"                 = "#e6550d",
  "TAM_like_MRC1"            = "#fd8d3c",
  "TAM_like_SPP1"            = "#fdae6b",
  "Monocytes_FCN1"           = "#fdd0a2",
  "Cycling_Myeloid"          = "#BCBD22",
  "Neutrophils"              = "#17BECF",
  "B_cells"                  = "#33a02c",
  "Plasma_cells"             = "#b2df8a",
  "CAF_like"                 = "#6a3d9a",
  "Endothelial_cells"        = "#b15928",
  "Mast_cells"               = "#fdbf6f",
  "Tumor_epithelial"         = "#a6cee3",
  "Tumor_epithelial_basal"   = "#8b0000",
  "Tumor_epithelial_AT2"     = "#c51b8a",
  "Tumor_epithelial_EMT"     = "#de2d26",
  "Ciliated_epithelial"      = "#74c476"
)

# Proportions
df <- seu_tme_mpr@meta.data %>%
  transmute(
    condition = PathResponse,
    cell_type = TME_cell_type
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
  scale_fill_manual(values = immuno_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  ylab("Cell-type proportion") +
  theme_classic() +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 13, face = "bold"),
    axis.text.y      = element_text(size = 11),
    axis.title.y     = element_text(size = 12),
    legend.position  = "right",
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(0.5, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.title     = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(fig_dir, "03_TME_Barplot_proportions_MPR_NMPR.png"),
       p_bar, width = 8, height = 6, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "03_TME_Barplot_proportions_MPR_NMPR.pdf"),
       p_bar, width = 8, height = 6, bg = "white")

cat("Done — TME barplot saved\n")