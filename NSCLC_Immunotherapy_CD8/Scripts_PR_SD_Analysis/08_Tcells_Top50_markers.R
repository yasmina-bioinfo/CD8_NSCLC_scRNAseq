#!/usr/bin/env Rscript
# ============================================================
# Immunotherapy_CD4_CD8_subset.R
# Subset CD4 et CD8 + top50 markers par cluster
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# -----------------------------
# Load
# -----------------------------
seu_t <- readRDS("Objects/05_Tcells_reclustered.rds")
DefaultAssay(seu_t) <- "RNA"

dir.create("Results/markers", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Définition des clusters
# -----------------------------
cd4_clusters <- c("1", "3", "4", "8", "11")
cd8_clusters <- c("0", "2", "5", "6", "7", "9", "10", "12", "13")

# -----------------------------
# Subset
# -----------------------------
seu_cd4 <- subset(seu_t, seurat_clusters %in% cd4_clusters)
seu_cd8 <- subset(seu_t, seurat_clusters %in% cd8_clusters)

cat("CD4 cells:", ncol(seu_cd4), "\n")
cat("CD8 cells:", ncol(seu_cd8), "\n")

# -----------------------------
# Fonction top50 markers
# -----------------------------
get_top50 <- function(seu_obj, label) {
  
  Idents(seu_obj) <- "seurat_clusters"
  
  markers <- FindAllMarkers(
    seu_obj,
    only.pos    = TRUE,
    min.pct     = 0.25,
    logfc.threshold = 0.25,
    test.use    = "wilcox"
  )
  
  top50 <- markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 50) %>%
    ungroup()
  
  write.csv(top50,
            paste0("Results/markers/top50_", label, ".csv"),
            row.names = FALSE)
  
  cat("\nTop50", label, "saved —", nrow(top50), "rows\n")
  return(top50)
}

# -----------------------------
# Calcul markers
# -----------------------------
cat("\nCalculating CD4 markers...\n")
top50_cd4 <- get_top50(seu_cd4, "CD4")

cat("\nCalculating CD8 markers...\n")
top50_cd8 <- get_top50(seu_cd8, "CD8")

# -----------------------------
# Sauvegarder les subsets
# -----------------------------
saveRDS(seu_cd4, "Objects/06_CD4_subset.rds")
saveRDS(seu_cd8, "Objects/06_CD8_subset.rds")

cat("\nDone — objects saved\n")