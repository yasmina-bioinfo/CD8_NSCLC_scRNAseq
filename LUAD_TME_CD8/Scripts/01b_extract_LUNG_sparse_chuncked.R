#!/usr/bin/env Rscript

# ============================================================
# GSE131907 — RAM-safe streaming conversion
# Read RAW (data.frame) -> keep only LUNG columns -> build dgCMatrix
# WITHOUT densifying (no as.matrix)
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
})

raw_path <- "Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
out_path <- "Data/GSE131907_Lung_only_UMI_matrix_sparse.rds"

message("[1] Starting script...")
message("[2] Reading RAW RDS...")
raw_mat <- readRDS(raw_path)
message("[3] RAW loaded. Class: ", paste(class(raw_mat), collapse=", "))
message("[4] RAW dim: ", paste(dim(raw_mat), collapse=" x "))

# Identify LUNG columns
message("[5] Finding LUNG columns...")
lung_j <- grep("LUNG", colnames(raw_mat))
message("[6] LUNG cells found: ", length(lung_j))

# Subset columns as data.frame (lightweight; still may take time)
message("[7] Subsetting RAW to LUNG columns...")
lung_df <- raw_mat[, lung_j, drop = FALSE]
rm(raw_mat); gc()
message("[8] LUNG df dim: ", paste(dim(lung_df), collapse=" x "))

# Build sparse dgCMatrix column-by-column (CSC format)
message("[9] Building sparse matrix (column-by-column)...")
n_genes <- nrow(lung_df)
n_cells <- ncol(lung_df)

i_list <- vector("list", n_cells)   # row indices per column
x_list <- vector("list", n_cells)   # values per column
p <- integer(n_cells + 1)           # column pointers (cumulative nnz)
p[1] <- 0L

# Progress every 500 columns (adjust if you want)
step <- 500L

for (j in seq_len(n_cells)) {
  v <- lung_df[[j]]
  
  # Make sure v is numeric (sometimes readRDS yields integer)
  # This avoids surprises in Matrix slots.
  if (!is.numeric(v)) v <- as.numeric(v)
  
  nz <- which(v != 0)
  i_list[[j]] <- nz - 1L           # 0-based for dgCMatrix
  x_list[[j]] <- v[nz]
  
  p[j + 1L] <- p[j] + length(nz)
  
  if (j %% step == 0L) {
    message("  ...processed columns: ", j, " / ", n_cells,
            " | nnz so far: ", p[j + 1L])
    gc()
  }
}

rm(lung_df); gc()

message("[10] Collapsing lists -> vectors (final assembly)...")
i <- as.integer(unlist(i_list, use.names = FALSE))
x <- as.numeric(unlist(x_list, use.names = FALSE))
rm(i_list, x_list); gc()

lung_sparse <- new("dgCMatrix",
                   i = i,
                   p = p,
                   x = x,
                   Dim = c(n_genes, n_cells),
                   Dimnames = list(NULL, NULL)
)

message("[11] Sparse built. Class: ", paste(class(lung_sparse), collapse=", "))
message("[12] Sparse dim: ", paste(dim(lung_sparse), collapse=" x "))
message("[13] nnz: ", length(lung_sparse@x))

saveRDS(lung_sparse, out_path)
message("[14] Saved -> ", out_path)