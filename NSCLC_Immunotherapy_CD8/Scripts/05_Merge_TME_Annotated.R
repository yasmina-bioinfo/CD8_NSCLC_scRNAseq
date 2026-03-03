#!/usr/bin/env Rscript

# ============================================================
# 04_TME_annotated.R
# Inject canonical cluster labels into TME object
# Input : 03_seu_umap_clustered.rds
# Output: 04_TME_annotated.rds
# Column: TME_cell_type
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

# -----------------------------
# Paths
# -----------------------------
in_rds  <- "objects/03_seu_umap_clustered.rds"
out_rds <- "objects/04_TME_annotated.rds"

# -----------------------------
# Load object
# -----------------------------
seu <- readRDS(in_rds)

# -----------------------------
# Canonical mapping (cluster -> label)
# -----------------------------
cluster_map <- c(
  "0"  = "T_Naive_Memory",
  "1"  = "Neutrophils",
  "2"  = "CD8_Cytotoxic_Exhausted",
  "3"  = "TAM_like",
  "4"  = "B_cells",
  "5"  = "TAM_like_MRC1",
  "6"  = "Tumor_epithelial",
  "7"  = "Plasma_cells",
  "8"  = "Monocytes_FCN1",
  "9"  = "Tumor_epithelial_basal",
  "10" = "Tumor_epithelial_AT2",
  "11" = "Tumor_epithelial_EMT",
  "12" = "Cycling_Tcells",
  "13" = "TAM_like_SPP1",
  "14" = "T_cells_Activated",
  "15" = "Mast_cells",
  "16" = "CAF_like",
  "17" = "Ciliated_epithelial",
  "18" = "Cycling_Myeloid",
  "19" = "Endothelial_cells"
)

# -----------------------------
# Safety checks: clusters exist
# -----------------------------
if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("ERROR: 'seurat_clusters' not found in meta.data. Found: ",
       paste(colnames(seu@meta.data), collapse = ", "))
}

cl <- as.character(seu$seurat_clusters)

missing_in_map <- setdiff(sort(unique(cl)), names(cluster_map))
if (length(missing_in_map) > 0) {
  stop("ERROR: These cluster IDs are missing in cluster_map: ",
       paste(missing_in_map, collapse = ", "))
}

# -----------------------------
# Inject annotation column
# -----------------------------
seu$TME_cell_type <- unname(cluster_map[cl])

# -----------------------------
# Validate: no NAs
# -----------------------------
na_n <- sum(is.na(seu$TME_cell_type))
if (na_n > 0) {
  stop("ERROR: TME_cell_type contains ", na_n, " NA values. Check cluster_map.")
}

message("\n=== Annotation injected: TME_cell_type ===")
print(table(seu$TME_cell_type))

# -----------------------------
# Save annotated object
# -----------------------------
saveRDS(seu, file = out_rds)
message("\nDONE: saved -> ", out_rds)