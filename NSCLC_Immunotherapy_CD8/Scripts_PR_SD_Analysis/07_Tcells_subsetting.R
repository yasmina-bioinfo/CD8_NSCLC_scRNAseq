#!/usr/bin/env Rscript

# ============================================================
# 05_Tcells_recluster.R
# Subset T cells from TME and recluster
# Input : 04_TME_annotated.rds
# Output: 05_Tcells_reclustered.rds + UMAP figures
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
in_rds  <- "objects/04_TME_annotated.rds"
out_rds <- "objects/05_Tcells_reclustered.rds"
figdir  <- "results/figures"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load TME object
# -----------------------------
seu <- readRDS(in_rds)

# -----------------------------
# Subset: all T cell types
# -----------------------------
t_labels <- c(
  "T_Naive_Memory",
  "CD8_Cytotoxic_Exhausted",
  "Cycling_Tcells",
  "T_cells_Activated"
)

stopifnot("TME_cell_type" %in% colnames(seu@meta.data))

seu_T <- subset(seu, subset = TME_cell_type %in% t_labels)

message("T cells subset dims: ", paste(dim(seu_T), collapse = " x "))
message("TME_cell_type composition:")
print(table(seu_T$TME_cell_type))

# -----------------------------
# Recluster cleanly
# -----------------------------
DefaultAssay(seu_T) <- "RNA"

# Optional but safe: clear old reductions/graphs
seu_T@reductions <- list()
seu_T@graphs <- list()

seu_T <- NormalizeData(seu_T, verbose = FALSE)
seu_T <- FindVariableFeatures(seu_T, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu_T <- ScaleData(seu_T, verbose = FALSE)
seu_T <- RunPCA(seu_T, npcs = 30, verbose = FALSE)

seu_T <- FindNeighbors(seu_T, dims = 1:30, verbose = FALSE)
seu_T <- FindClusters(seu_T, resolution = 0.5, verbose = FALSE)
seu_T <- RunUMAP(seu_T, dims = 1:30, verbose = FALSE)

# -----------------------------
# Quick UMAP (by clusters)
# -----------------------------
p_umap_T <- DimPlot(seu_T, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.2) +
  ggtitle("T cells reclustering") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(figdir, "Tcells_UMAP_reclustered.png"), p_umap_T,
       width = 8, height = 6, dpi = 450, bg = "white")
ggsave(file.path(figdir, "Tcells_UMAP_reclustered.pdf"), p_umap_T,
       width = 8, height = 6, bg = "white")

# -----------------------------
# Save object
# -----------------------------
saveRDS(seu_T, out_rds)
message("DONE: saved -> ", out_rds)
message("DONE: figures -> ", figdir)