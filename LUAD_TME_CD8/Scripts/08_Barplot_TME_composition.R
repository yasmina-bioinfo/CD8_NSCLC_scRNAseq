# 08_barplot_TME_composition.R
# Stacked barplot of major cell types per Sample_Origin (proportions)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- Load object ----
obj <- readRDS("objects/07_seurat_with_sample_origin.rds")

# ---- Palette (same as UMAP) ----
tme_colors <- c(
  "CD8 T cells"                = "#1f78b4",
  "CD4 T cells"                = "#6baed6",
  "B lineage"                  = "#33a02c",
  "Macrophages"                = "#ff7f00",
  "Monocytes"                  = "#fdae6b",
  "Dendritic cells"            = "#fb9a99",
  "Mast cells"                 = "#fdbf6f",
  "Fibroblasts"                = "#6a3d9a",
  "Endothelial cells"          = "#b15928",
  "Tumor basal/proliferative"  = "#67000d",
  "Tumor hypoxic/remodeling"   = "#cb181d",
  "Tumor secretory/mucinous"   = "#f03b20",
  "Lung epithelial"            = "#9ecae1",
  "Other"                      = "#bdbdbd"
)

# ---- Prepare data ----
# Table counts
df <- obj@meta.data %>%
  dplyr::count(Sample_Origin, major_celltype) %>%
  group_by(Sample_Origin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Biological order 
compartment_order <- c("tLung", "nLung", "PE", "tL/B", "nLN", "mLN", "mBrain")
df$Sample_Origin <- factor(df$Sample_Origin, levels = compartment_order)

# ---- Barplot ----
p_bar <- ggplot(df, aes(x = Sample_Origin, y = prop, fill = major_celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = tme_colors) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "TME composition across compartments",
    x = "Compartment",
    y = "Proportion of cells"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

  guides(fill = guide_legend(ncol = 2))

# ---- Save ----
ggsave("results/figures/Barplot_TME_composition.png",
       p_bar, width = 9, height = 6, dpi = 400, bg = "white")

ggsave("results/figures/Barplot_TME_composition.pdf",
       p_bar, width = 9, height = 6, bg = "white")


message("Barplot saved.")