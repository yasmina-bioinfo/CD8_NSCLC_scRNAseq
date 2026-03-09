#!/usr/bin/env Rscript

# ============================================================
# 09_Subset_T_Cells.R
# Subset T cells -> subset nLung/tLung -> save object
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

# -----------------------------
# Paths (adjust if needed)
# -----------------------------
in_rds  <- "objects/06_seurat_TME_with_SampleOrigin_and_orderedCelltypes.rds"                 # <- objet LUAD post-QC/PCA
out_rds <- "objects/09_Tcells_nLung_tLung_preReclustered.rds" # <- output

# -----------------------------
# Load
# -----------------------------
seu <- readRDS(in_rds)

# -----------------------------
# Sanity checks
# -----------------------------
if (!"Cell_type" %in% colnames(seu@meta.data)) {
  stop("Meta.data column 'Cell_type' not found. Check your object metadata.")
}
if (!"Sample_Origin" %in% colnames(seu@meta.data)) {
  stop("Meta.data column 'Sample_Origin' not found. Check your object metadata.")
}

# -----------------------------
# Subset: T cells
# -----------------------------
Idents(seu) <- "major_celltype"
seu_T <- subset(seu, idents = c("T lymphocytes", "T/NK cells"))

# -----------------------------
# Subset: keep only nLung / tLung
# -----------------------------
seu_T <- subset(seu_T, subset = Sample_Origin %in% c("nLung", "tLung"))

# -----------------------------
# Minimal outputs for validation
# -----------------------------
cat("\n=== 09_Subset_T_Cells.R : VALIDATION ===\n")
cat("Dims (features x cells): ", paste(dim(seu_T), collapse = " x "), "\n\n")
cat("Sample_Origin counts:\n")
print(table(seu_T$Sample_Origin))

# -----------------------------
# Save
# -----------------------------
dir.create("Objects", showWarnings = FALSE, recursive = TRUE)
saveRDS(seu_T, out_rds)
cat("\nSaved -> ", out_rds, "\n")