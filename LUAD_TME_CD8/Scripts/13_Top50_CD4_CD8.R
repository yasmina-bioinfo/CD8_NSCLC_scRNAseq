#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

dir.create("Results/markers", showWarnings = FALSE, recursive = TRUE)
dir.create("Objects", showWarnings = FALSE, recursive = TRUE)

recluster_and_top50 <- function(in_rds, prefix, dims_use = 1:15, resolution = 0.5) {
  
  seu <- readRDS(in_rds)
  DefaultAssay(seu) <- "RNA"
  
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = 2000)
  seu <- ScaleData(seu, features = VariableFeatures(seu))
  seu <- RunPCA(seu, npcs = 30, features = VariableFeatures(seu))
  
  seu <- FindNeighbors(seu, dims = dims_use)
  seu <- FindClusters(seu, resolution = resolution)
  
  Idents(seu) <- "seurat_clusters"
  
  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  top50 <- markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 50) %>%
    ungroup()
  
  write.table(
    top50,
    file = file.path("Results/markers", paste0("14_", prefix, "_Top50Markers.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  out_rds <- file.path("Objects", paste0("14_", prefix, "_reclustered.rds"))
  saveRDS(seu, out_rds)
  
  cat("\n== ", prefix, " ==\n", sep = "")
  cat("Cells: ", ncol(seu), "\n", sep = "")
  cat("Clusters:\n")
  print(table(seu$seurat_clusters))
  cat("Saved object -> ", out_rds, "\n", sep = "")
  cat("Saved Top50 -> ", file.path("Results/markers", paste0("14_", prefix, "_Top50Markers.csv")), "\n", sep = "")
}

recluster_and_top50("Objects/13_CD4_raw_subset.rds", "CD4", dims_use = 1:15, resolution = 0.5)
recluster_and_top50("Objects/13_CD8_raw_subset.rds", "CD8", dims_use = 1:15, resolution = 0.5)