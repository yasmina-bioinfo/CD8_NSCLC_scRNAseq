# 07_major_celltype_mapping.R
# Purpose: Create a reproducible "major_celltype" column from Seurat cluster identities (res=0.5)

suppressPackageStartupMessages({
  library(Seurat)
})

# -------------------------
# Paths (adapt to your repo)
# -------------------------
in_rds  <- "objects/05_seurat_with_umap_clusters.rds"
out_rds <- "objects/06_seurat_with_major_celltype.rds"
out_dir <- "Results/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load object
# -------------------------
obj <- readRDS(in_rds)

# -------------------------
# Ensure correct identities
# -------------------------
# If your clusters are stored under another name, change it here.
# Common options: "seurat_clusters" or "RNA_snn_res.0.5"
if ("RNA_snn_res.0.5" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "RNA_snn_res.0.5"
} else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "seurat_clusters"
} else {
  stop("No clustering identity found. Expected 'RNA_snn_res.0.5' or 'seurat_clusters' in metadata.")
}

# Convert to numeric-like levels if needed
clusters <- levels(Idents(obj))
message("Clusters detected: ", paste(clusters, collapse = ", "))

# -------------------------
# Create major_celltype
# -------------------------
obj$major_celltype <- "Unassigned"

# CD4 T cells
obj$major_celltype[Idents(obj) %in% c("0","9","29")] <- "CD4 T cells"

# CD8 T cells
obj$major_celltype[Idents(obj) %in% c("1","3","21")] <- "CD8 T cells"

# B lineage (B + Plasma)
obj$major_celltype[Idents(obj) %in% c("2","16")] <- "B lineage"

# Myeloid
obj$major_celltype[Idents(obj) %in% c("6")] <- "Monocytes"
obj$major_celltype[Idents(obj) %in% c("4","5","27")] <- "Macrophages"
obj$major_celltype[Idents(obj) %in% c("14","18")] <- "Dendritic cells"
obj$major_celltype[Idents(obj) %in% c("11")] <- "Mast cells"

# Stromal
obj$major_celltype[Idents(obj) %in% c("8")] <- "Fibroblasts"
obj$major_celltype[Idents(obj) %in% c("17")] <- "Endothelial cells"

# Tumor epithelial (3 strategic groups)
obj$major_celltype[Idents(obj) %in% c("13","20","23")] <- "Tumor basal/proliferative"
obj$major_celltype[Idents(obj) %in% c("12","22")] <- "Tumor hypoxic/remodeling"
obj$major_celltype[Idents(obj) %in% c("7","19","24","26")] <- "Tumor secretory/mucinous"

# Other epithelial
obj$major_celltype[Idents(obj) %in% c("10","28")] <- "Lung epithelial"

# Other (Cycling + Rare)
obj$major_celltype[Idents(obj) %in% c("15","25")] <- "Other"

# -------------------------
# QC checks (must pass)
# -------------------------
tab_major <- sort(table(obj$major_celltype), decreasing = TRUE)
print(tab_major)

if (any(obj$major_celltype == "Unassigned")) {
  unassigned_clusters <- sort(unique(as.character(Idents(obj)[obj$major_celltype == "Unassigned"])))
  stop(
    "Some cells remain Unassigned. Missing cluster(s): ",
    paste(unassigned_clusters, collapse = ", "),
    "\nUpdate the mapping above to include them."
  )
}

# -------------------------
# Export tables for Excel
# -------------------------
write.csv(
  as.data.frame(tab_major),
  file = file.path(out_dir, "major_celltype_counts.csv"),
  row.names = FALSE
)

# Optional: also export cluster -> major mapping as a quick reference
cluster_major_map <- data.frame(
  cluster = levels(Idents(obj)),
  major_celltype = sapply(levels(Idents(obj)), function(cl){
    # take the most frequent label in that cluster
    names(sort(table(obj$major_celltype[Idents(obj) == cl]), decreasing = TRUE))[1]
  })
)
write.csv(
  cluster_major_map,
  file = file.path(out_dir, "cluster_to_major_celltype.csv"),
  row.names = FALSE
)

# -------------------------
# Save updated object
# -------------------------
saveRDS(obj, out_rds)
message("Saved: ", out_rds)
message("Exported tables to: ", out_dir)