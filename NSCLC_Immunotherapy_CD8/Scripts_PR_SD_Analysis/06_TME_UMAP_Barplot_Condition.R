#!/usr/bin/env Rscript

# ============================================================
# 05_TME_visualize.R
# UMAP (TME_cell_type) + Barplot split PR/SD
# Input : 04_TME_annotated.rds
# Output: Results/figures/*.png + *.pdf
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# Paths
# -----------------------------
in_rds <- "Objects/04_TME_annotated.rds"
outdir <- "results/figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load object
# -----------------------------
seu <- readRDS(in_rds)

# -----------------------------
# REQUIRED metadata columns
# -----------------------------
celltype_col <- "TME_cell_type"

# >>> CHANGE THIS if your column name differs <<<
cond_col <- "RECIST"   # expected values: PR / SD

stopifnot(celltype_col %in% colnames(seu@meta.data))
stopifnot(cond_col %in% colnames(seu@meta.data))

# -----------------------------
# Colors (your palette)
# -----------------------------
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

# Ensure consistent legend order + avoid missing colors
seu[[celltype_col]] <- factor(seu[[celltype_col]][, 1], levels = names(immuno_colors))

# -----------------------------
# UMAP - overall
# -----------------------------
p_umap <- DimPlot(
  seu,
  reduction = "umap",
  group.by = celltype_col,
  cols = immuno_colors,
  label = FALSE,
  pt.size = 0.15
) +
  ggtitle("TME cell types") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.7, "cm"),
    legend.spacing.y = unit(0.4, "cm")
  )

ggsave(file.path(outdir, "TME_UMAP_TME_cell_type.png"), p_umap,
       width = 8, height = 8, dpi = 450, bg = "white")
ggsave(file.path(outdir, "TME_UMAP_TME_cell_type.pdf"), p_umap,
       width = 8, height = 8, bg = "white")

# -----------------------------
# UMAP - split by PR/SD
# -----------------------------
p_umap_split <- DimPlot(
  seu,
  reduction = "umap",
  group.by = celltype_col,
  split.by = cond_col,
  cols = immuno_colors,
  label = FALSE,
  pt.size = 0.15
) +
  ggtitle("TME cell types") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 9)
  )

ggsave(file.path(outdir, "TME_UMAP_TME_cell_type_split_PR_SD.png"), p_umap_split,
       width = 10, height = 5, dpi = 450, bg = "white")
ggsave(file.path(outdir, "TME_UMAP_TME_cell_type_split_PR_SD.pdf"), p_umap_split,
       width = 10, height = 5, bg = "white")

# -----------------------------
# Barplot: proportions by condition (PR/SD)
# -----------------------------
df <- seu@meta.data %>%
  transmute(
    condition = .data[[cond_col]],
    cell_type = .data[[celltype_col]]
  ) %>%
  filter(!is.na(condition), !is.na(cell_type))

df_prop <- df %>%
  count(condition, cell_type, name = "n") %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar <- ggplot(df_prop, aes(x = condition, y = prop, fill = cell_type)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = immuno_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Cell-type proportion") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 13, face = "bold"),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    
    legend.position = "right",
    legend.text  = element_text(size = 11),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.title = element_blank()
  
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(outdir, "TME_BARPLOT_proportions_PR_SD.png"), p_bar,
       width = 8, height = 6, dpi = 450, bg = "white")
ggsave(file.path(outdir, "TME_BARPLOT_proportions_PR_SD.pdf"), p_bar,
       width = 8, height = 6, bg = "white")

message("DONE: UMAP + barplot saved in -> ", outdir)