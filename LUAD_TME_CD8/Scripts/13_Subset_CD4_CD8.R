#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

seu <- readRDS("Objects/12_Tcells_clean_annotated.rds")
DefaultAssay(seu) <- "RNA"

# Retirer NK_like si tu veux rester strictement T classique
seu <- subset(seu, subset = Tcell_state != "NK_like_Cytotoxic")

# Expression matrix
expr <- GetAssayData(seu, layer = "data")

cd4_pos <- expr["CD4", ] > 0
cd8_pos <- expr["CD8A", ] > 0 | expr["CD8B", ] > 0

cd4_cells <- colnames(seu)[cd4_pos & !cd8_pos]
cd8_cells <- colnames(seu)[cd8_pos]

seu_cd4 <- subset(seu, cells = cd4_cells)
seu_cd8 <- subset(seu, cells = cd8_cells)

cat("\nCD4 cells:", ncol(seu_cd4))
cat("\nCD8 cells:", ncol(seu_cd8), "\n")

saveRDS(seu_cd4, "Objects/13_CD4_raw_subset.rds")
saveRDS(seu_cd8, "Objects/13_CD8_raw_subset.rds")