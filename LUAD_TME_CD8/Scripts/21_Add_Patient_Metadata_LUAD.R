#!/usr/bin/env Rscript
# ============================================================
# 21_Add_Patient_Metadata_LUAD.R
# Recover patient IDs from cell barcodes (lost during downsampling)
# Transfer to CD8 annotated object
# ============================================================
library(Seurat)

object_dir <- "Objects"

# -----------------------------
# Load raw object (patient info encoded in cell names)
# -----------------------------
seu_raw <- readRDS("Objects/01_seurat_40k_raw_sparse.rds")

# -----------------------------
# Extract patient ID from cell barcodes
# Pattern: BARCODE_CONDITION_PATIENTID
# e.g. GGCGTGTCACGGTAAG_EFFUSION_13 → EFFUSION_13
# -----------------------------
patient_map <- sub(".*_([^_]+_[^_]+)$", "\\1", colnames(seu_raw))
names(patient_map) <- colnames(seu_raw)

cat("Patient IDs found:\n")
print(table(patient_map))

# -----------------------------
# Load CD8 annotated object
# -----------------------------
seu_luad <- readRDS("Objects/15_CD8_annotated_states.rds")

# -----------------------------
# Transfer patient metadata
# -----------------------------
matching_cells <- intersect(colnames(seu_luad), names(patient_map))
cat("\nMatching cells:", length(matching_cells), "/", ncol(seu_luad), "\n")

seu_luad$Patient <- NA
seu_luad$Patient[matching_cells] <- patient_map[matching_cells]

cat("\nPatient distribution in CD8 object:\n")
print(table(seu_luad$Patient, useNA = "always"))

# -----------------------------
# Save updated object
# -----------------------------
saveRDS(seu_luad, "Objects/16_CD8_LUAD_with_Patient.rds")

cat("\nDone — 16_CD8_LUAD_with_Patient.rds saved\n")
