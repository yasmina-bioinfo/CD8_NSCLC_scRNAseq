#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Script 02: QC (violin + scatter) + filtering
# Input : objects/01_seu_raw.rds
# Output: objects/02_seu_qc.rds + results/figures/QC_*.png
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# Paths (EDIT THIS)
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy"
IN_OBJ   <- file.path(DATA_DIR, "objects", "01_seu_raw.rds")
OUT_OBJ  <- file.path(DATA_DIR, "objects")
OUT_FIG  <- file.path(DATA_DIR, "results", "figures")

dir.create(OUT_OBJ, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load raw object
# -----------------------------
seu <- readRDS(IN_OBJ)
DefaultAssay(seu) <- "RNA"

# -----------------------------
# 2) QC metrics
# -----------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# -----------------------------
# 3) Violin QC (by RECIST)
# -----------------------------
p_vln <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "RECIST",
  ncol = 3,
  pt.size = 0
) & theme_bw()

ggsave(
  filename = file.path(OUT_FIG, "QC_violin_by_RECIST.png"),
  plot = p_vln,
  width = 12, height = 4, dpi = 300, bg = "white"
)

# -----------------------------
# 4) Scatter QC (classic)
# -----------------------------
p_sc1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme_bw()
p_sc2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_bw()
p_sc  <- p_sc1 + p_sc2

ggsave(
  filename = file.path(OUT_FIG, "QC_scatter.png"),
  plot = p_sc,
  width = 10, height = 4, dpi = 300, bg = "white"
)

# -----------------------------
# 5) Filtering (conservative defaults)
#    -> You can adjust after inspecting QC plots
# -----------------------------
min_features <- 300
max_features <- 6500
max_mt <- 18

seu_qc <- subset(seu, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_mt)

message("After QC filtering:")
print(seu_qc)
message("RECIST distribution after QC:")
print(table(seu_qc$RECIST, useNA = "ifany"))

# -----------------------------
# 6) Save QC object
# -----------------------------
saveRDS(seu_qc, file.path(OUT_OBJ, "02_seu_qc.rds"))
message("Saved: objects/02_seu_qc.rds")
message("DONE Script 02")