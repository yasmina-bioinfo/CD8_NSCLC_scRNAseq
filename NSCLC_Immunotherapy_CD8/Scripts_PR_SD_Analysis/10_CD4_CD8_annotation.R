# ============================================================
# 07_CD4_CD8_annotation.R
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

seu_cd4 <- readRDS("objects/06_CD4_subset.rds")
seu_cd8 <- readRDS("objects/06_CD8_subset.rds")

# -----------------------------
# Mapping annotations
# -----------------------------
cd4_map <- c(
  "1"  = "CD4_Naive_CentralMemory",
  "3"  = "CD4_Tfh_like",
  "4"  = "Treg_Activated",
  "8"  = "CD4_Early_Activated_NR4A_high",
  "11" = "Treg_Effector"
)

cd8_map <- c(
  "0"  = "CD8_Effector_GZMK",
  "2"  = "CD8_Exhausted_Terminal",
  "5"  = "CD8_Terminal_CX3CR1",
  "6"  = "CD8_Proliferating",
  "7"  = "CD8_TRM_like",
  "9"  = "CD8_Early_Activated_NR4A_high",
  "10" = "CD8_Activated_HLAII_high",
  "12" = "CD8_IFN_Stress_Response"
)

# -----------------------------
# Appliquer le mapping
# -----------------------------
seu_cd4$cell_state <- dplyr::recode(
  as.character(seu_cd4$seurat_clusters),
  !!!cd4_map
)

seu_cd8$cell_state <- dplyr::recode(
  as.character(seu_cd8$seurat_clusters),
  !!!cd8_map
)

# Définir comme identité active
Idents(seu_cd4) <- "cell_type"
Idents(seu_cd8) <- "cell_type"

# -----------------------------
# Vérification
# -----------------------------
cat("=== CD4 annotation ===\n")
print(table(seu_cd4$cell_state))

cat("\n=== CD8 annotation ===\n")
print(table(seu_cd8$cell_state))

cat("\n=== CD4 response distribution ===\n")
print(table(seu_cd4$cell_state, seu_cd4$response))

cat("\n=== CD8 response distribution ===\n")
print(table(seu_cd8$cell_state, seu_cd8$response))

# -----------------------------
# Sauvegarder
# -----------------------------
saveRDS(seu_cd4, "objects/07_CD4_annotated.rds")
saveRDS(seu_cd8, "objects/07_CD8_annotated.rds")

cat("\nDone — objects saved\n")