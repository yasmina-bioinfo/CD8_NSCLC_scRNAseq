#!/usr/bin/env Rscript
# ============================================================
# 07_CellChat_MPR_NMPR_Create.R
# Create CellChat objects — MPR vs NMPR
# Anti-PD1 neoadjuvant context — focus on CD8 checkpoint interactions
# ============================================================
suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(dplyr)
})

fig_dir    <- "Results/figures"
object_dir <- "Objects"
dir.create(fig_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(object_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load full TME object (with PathResponse metadata)
# -----------------------------
seu <- readRDS("Objects/04_TME_MPR_NMPR.rds")
DefaultAssay(seu) <- "RNA"

# -----------------------------
# Helper — create one CellChat object
# -----------------------------
create_cc <- function(seu_obj, condition_label) {
  
  cat("Creating CellChat for:", condition_label, "\n")
  
  data_input <- GetAssayData(seu_obj, assay = "RNA", layer = "data")
  meta       <- seu_obj@meta.data
  
  cc <- createCellChat(
    object   = data_input,
    meta     = meta,
    group.by = "TME_cell_type"
  )
  
  # Secreted Signaling + Cell-Cell Contact
  # captures PD1-PDL1, TIGIT, TIM3 axes
  CellChatDB.use <- subsetDB(
    CellChatDB.human,
    search = c("Secreted Signaling", "Cell-Cell Contact")
  )
  cc@DB <- CellChatDB.use
  
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  
  cc <- computeCommunProb(cc, type = "triMean", nboot = 20)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  
  return(cc)
}

# -----------------------------
# Subset MPR and NMPR
# -----------------------------
seu_MPR  <- subset(seu, PathResponse == "MPR")
seu_NMPR <- subset(seu, PathResponse == "NMPR")

# Add samples column (required by CellChat for multi-sample awareness)
seu_MPR@meta.data$samples  <- seu_MPR@meta.data$Patient
seu_NMPR@meta.data$samples <- seu_NMPR@meta.data$Patient

cat("MPR cells:",  ncol(seu_MPR),  "\n")
cat("NMPR cells:", ncol(seu_NMPR), "\n")

# -----------------------------
# Create objects
# -----------------------------
cc_MPR  <- create_cc(seu_MPR,  "MPR")
cc_NMPR <- create_cc(seu_NMPR, "NMPR")

# -----------------------------
# Merge
# -----------------------------
cc_merged <- mergeCellChat(
  list(cc_MPR, cc_NMPR),
  add.names = c("MPR", "NMPR")
)

# -----------------------------
# Save
# -----------------------------
saveRDS(cc_MPR,    "Objects/07_CellChat_MPR.rds")
saveRDS(cc_NMPR,   "Objects/07_CellChat_NMPR.rds")
saveRDS(cc_merged, "Objects/07_CellChat_merged.rds")

cat("\nDone — CellChat MPR/NMPR objects saved\n")
