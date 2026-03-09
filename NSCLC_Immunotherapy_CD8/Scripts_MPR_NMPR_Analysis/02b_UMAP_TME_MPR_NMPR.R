# ============================================================
# 02b_UMAP_TME_MPR_NMPR.R
# UMAP TME colored by cell type and by condition (MPR/NMPR)
# ============================================================

library(Seurat)
library(ggplot2)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

seu_tme_mpr <- readRDS("Objects/04_TME_MPR_NMPR.rds")

# Colors
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

# UMAP by cell type
p1 <- DimPlot(seu_tme_mpr, group.by = "TME_cell_type",
              cols = immuno_colors, pt.size = 0.3,
              label = FALSE) +
  ggtitle("TME cell types — MPR/NMPR") +
  theme(plot.title = element_text(size = 13, face = "bold"))

# UMAP by condition
p2 <- DimPlot(seu_tme_mpr, group.by = "PathResponse",
              cols = c("MPR" = "#4575B4", "NMPR" = "#D73027"),
              pt.size = 0.3) +
  ggtitle("Condition — MPR vs NMPR") +
  theme(plot.title = element_text(size = 13, face = "bold"))

# Save
ggsave(file.path(fig_dir, "02b_UMAP_TME_celltypes_MPR_NMPR.png"),
       p1, width = 10, height = 7, dpi = 450, bg = "white")

ggsave(file.path(fig_dir, "02b_UMAP_TME_condition_MPR_NMPR.png"),
       p2, width = 8, height = 6, dpi = 450, bg = "white")

cat("Done — TME UMAPs saved\n")