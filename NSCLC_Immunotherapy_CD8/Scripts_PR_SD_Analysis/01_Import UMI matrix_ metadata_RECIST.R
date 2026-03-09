#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Script 01: Import UMI matrix + metadata (RECIST)
# Output: objects/01_seu_raw.rds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(readxl)
  library(Matrix)
})

# -----------------------------
# Paths (EDIT THIS)
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy"
OUT_OBJ  <- file.path(DATA_DIR, "objects")
dir.create(OUT_OBJ, showWarnings = FALSE, recursive = TRUE)

UMI_FILE  <- file.path(DATA_DIR, "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
META_FILE <- file.path(DATA_DIR, "GSE207422_NSCLC_scRNAseq_metadata.xlsx")

# -----------------------------
# 1) Read UMI matrix (genes x cells)
# -----------------------------
message("Reading UMI matrix...")
counts_dt <- fread(UMI_FILE)

stopifnot("Gene" %in% colnames(counts_dt))
gene_names <- counts_dt$Gene
counts_dt$Gene <- NULL

counts_mat <- as.matrix(counts_dt)
rownames(counts_mat) <- gene_names

rm(counts_dt); gc()

message("Counts dim (genes x cells): ", paste(dim(counts_mat), collapse = " x "))

# -----------------------------
# 2) Create Seurat object (raw)
# -----------------------------
seu <- CreateSeuratObject(counts = counts_mat, project = "GSE207422", min.cells = 3, min.features = 200)

# -----------------------------
# 3) Read sample-level metadata + add RECIST per cell
# -----------------------------
meta <- read_excel(META_FILE)

required_cols <- c("Sample", "RECIST")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) stop("Metadata missing: ", paste(missing_cols, collapse = ", "))

# Cell names: BD_immune01_612637  -> Sample: BD_immune01
seu$Sample <- sub("(_[0-9]+)$", "", colnames(seu))

if (!all(seu$Sample %in% meta$Sample)) {
  bad <- unique(seu$Sample[!(seu$Sample %in% meta$Sample)])
  stop("Some Samples not found in metadata. Example: ", paste(head(bad, 5), collapse = ", "))
}

recist_map <- setNames(meta$RECIST, meta$Sample)
seu$RECIST <- recist_map[seu$Sample]

message("RECIST distribution (cells):")
print(table(seu$RECIST, useNA = "ifany"))

# -----------------------------
# 4) Save raw object
# -----------------------------
saveRDS(seu, file.path(OUT_OBJ, "01_seu_raw.rds"))
message("Saved: objects/01_seu_raw.rds")
message("DONE Script 01")