library(Seurat)

# ---- Load 40k matrix (Matrix Market format) ----
counts_sparse <- ReadMtx(
  mtx = "Data/subsample_40k/umi_counts_40k.mtx",
  features = "Data/subsample_40k/features_40k.tsv",
  cells = "Data/subsample_40k/barcodes_40k.tsv",
  feature.column = 1
)

# ---- Create Seurat object ----
seu <- CreateSeuratObject(
  counts = counts_sparse,
  project = "LUAD_40k"
)

# ---- Save object ----
dir.create("objects", showWarnings = FALSE)
saveRDS(seu, "objects/01_seurat_40k_raw_sparse.rds")

# ---- Quick validation ----
print(dim(seu))
print(GetAssayData(seu, assay = "RNA", layer = "counts")[1:3, 1:3])

