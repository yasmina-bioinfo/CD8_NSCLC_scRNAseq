# 07_add_sample_origin_from_GEO.R
# Add Sample_Origin (tLung/nLung/mLN/nLN/mBrain/PE/tL-B...) from GEO annotation file
# Robust approach: write directly into obj@meta.data using exact cell names

suppressPackageStartupMessages({
  library(Seurat)
})

# ---- Paths (adapt if needed) ----
seurat_in  <- "objects/06_seurat_with_major_celltype.rds"
geo_anno   <- "data/GSE131907_Lung_Cancer_cell_annotation.txt"
seurat_out <- "objects/07_seurat_with_sample_origin.rds"

# ---- Load ----
obj  <- readRDS(seurat_in)
anno <- read.delim(geo_anno, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- Build mapping: core barcode -> Sample_Origin ----
origin_map <- setNames(anno$Sample_Origin, anno$Barcode)

# ---- Extract core barcode from Seurat cell names ----
# Example: "GGCGTGTCACGGTAAG_EFFUSION_13" -> "GGCGTGTCACGGTAAG"
seu_core <- sub("_.*$", "", colnames(obj))

# ---- Create vector aligned to Seurat cell names ----
v <- origin_map[seu_core]
names(v) <- colnames(obj)

# ---- Inject into metadata (guaranteed alignment) ----
obj@meta.data$Sample_Origin <- v[rownames(obj@meta.data)]

# ---- QC checks ----
stopifnot(nrow(obj@meta.data) == ncol(obj))
if (any(is.na(obj@meta.data$Sample_Origin))) {
  na_n <- sum(is.na(obj@meta.data$Sample_Origin))
  stop("Sample_Origin mapping produced NA for ", na_n, " cells. Check barcode formats.")
}

message("\nSample_Origin counts:")
print(table(obj@meta.data$Sample_Origin))

# Optional: quick cross-check with major cell types (should not error)
if ("major_celltype" %in% colnames(obj@meta.data)) {
  message("\nCross-tab Sample_Origin x major_celltype (head):")
  xt <- table(obj@meta.data$Sample_Origin, obj@meta.data$major_celltype)
  print(xt[1:min(5,nrow(xt)), 1:min(6,ncol(xt))])
}

# ---- Save ----
saveRDS(obj, "objects/07_seurat_with_sample_origin.rds")
saveRDS(obj, seurat_out)
message("\nSaved: ", seurat_out)