# 07_umap_TME_plots.R
# Generate UMAPs for TME with consistent colors
# - UMAP colored by major_celltype
# - UMAP colored by Cell_type.refined (optional)
# - UMAP split by Sample_Origin (major_celltype colors)

suppressPackageStartupMessages({
  library(Seurat)
})

# ---- Paths ----
seurat_in <- "objects/07_seurat_with_sample_origin.rds"
fig_dir   <- "Results/figures"

# ---- Load ----
obj <- readRDS(seurat_in)

# ---- Palette (Tumor vs Immune) ----
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

# ---- Safety checks ----
stopifnot("major_celltype" %in% colnames(obj@meta.data))
missing_cols <- setdiff(unique(obj$major_celltype), names(tme_colors))
if (length(missing_cols) > 0) {
  stop("Palette missing these major_celltype labels: ", paste(missing_cols, collapse = ", "))
}
stopifnot("Sample_Origin" %in% colnames(obj@meta.data))

# Ensure output folder exists
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---- Ensure Sample_Origin is populated (rebuild if NA) ----
if ("Sample_Origin" %in% colnames(obj@meta.data) && all(is.na(obj@meta.data$Sample_Origin))) {
  
  anno <- read.delim("data/GSE131907_Lung_Cancer_cell_annotation.txt",
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  origin_map <- setNames(anno$Sample_Origin, anno$Barcode)
  seu_core <- sub("_.*$", "", colnames(obj))
  
  v <- origin_map[seu_core]
  names(v) <- colnames(obj)
  
  obj@meta.data$Sample_Origin <- v[rownames(obj@meta.data)]
  
  stopifnot(sum(is.na(obj@meta.data$Sample_Origin)) == 0)
  
  # Save fixed object so it persists
  saveRDS(obj, "objects/07_seurat_with_sample_origin.rds")
  message("Fixed and saved Sample_Origin into objects/07_seurat_with_sample_origin.rds")
}


# ---- 1) UMAP: major_celltype (no labels, legend only) ----
p_major <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "major_celltype",
  cols = tme_colors,
  label = FALSE
) +
  ggplot2::ggtitle("Tumor Microenvironment — Major Cell Types") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 11)
  )

ggplot2::ggsave(
  filename = file.path(fig_dir, "UMAP_major_celltype.png"),
  plot = p_major,
  width = 8,
  height = 6,
  dpi = 400,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(fig_dir, "UMAP_major_celltype.pdf"),
  plot = p_major,
  width = 8,
  height = 6,
  bg = "white"
)

# ---- 2) UMAP: Cell_type.refined (if present) ----
if ("Cell_type.refined" %in% colnames(obj@meta.data)) {
  p_ref <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "Cell_type.refined",
    label = FALSE
  ) + ggplot2::ggtitle("TME UMAP — refined cell types")
  
  ggplot2::ggsave(
    filename = file.path(fig_dir, "UMAP_celltype_refined.png"),
    plot = p_ref, width = 9, height = 6, dpi = 400, bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(fig_dir, "UMAP_celltype_refined.pdf"),
    plot = p_ref, width = 9, height = 6, bg = "white"
  )
}

# ---- 3) UMAP split by Sample_Origin (horizontal grid) ----
p_split <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "major_celltype",
  cols = tme_colors,
  split.by = "Sample_Origin",
  label = FALSE,
  ncol = length(unique(obj$Sample_Origin))
) +
  ggplot2::ggtitle("Major Cell Types across Compartments (Sample_Origin)") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 10)
  )

ggplot2::ggsave(
  filename = file.path(fig_dir, "UMAP_major_celltype_split_Sample_Origin.png"),
  plot = p_split, width = 16, height = 6, dpi = 400, bg = "white"
)

ggplot2::ggsave(
  filename = file.path(fig_dir, "UMAP_major_celltype_split_Sample_Origin.pdf"),
  plot = p_split, width = 16, height = 6, bg = "white"
)


message("Done. Figures saved to: ", fig_dir)