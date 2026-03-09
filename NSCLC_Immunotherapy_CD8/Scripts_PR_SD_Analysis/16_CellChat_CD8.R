#!/usr/bin/env Rscript
# ============================================================
# 12_CD8_CellChat_create.R
# Create CellChat objects — PR vs SD
# Anti-PD1 context — focus on CD8 checkpoint interactions
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
# Load full TME object
# -----------------------------
seu <- readRDS("Objects/04_TME_annotated.rds")
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
  
  # nboot = 0 to reduce computation time
  cc <- computeCommunProb(cc, type = "triMean", nboot = 20)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  
  return(cc)
}

# -----------------------------
# Subset PR and SD
# -----------------------------
seu_PR <- subset(seu, RECIST == "PR")
seu_SD <- subset(seu, RECIST == "SD")

cat("PR cells:", ncol(seu_PR), "\n")
cat("SD cells:", ncol(seu_SD), "\n")

# -----------------------------
# Create objects
# -----------------------------
cc_PR <- create_cc(seu_PR, "PR")
cc_SD <- create_cc(seu_SD, "SD")

# -----------------------------
# Merge
# -----------------------------
cc_merged <- mergeCellChat(
  list(cc_PR, cc_SD),
  add.names = c("PR", "SD")
)

# -----------------------------
# Save
# -----------------------------
saveRDS(cc_PR,     "Objects/12_CellChat_PR.rds")
saveRDS(cc_SD,     "Objects/12_CellChat_SD.rds")
saveRDS(cc_merged, "Objects/12_CellChat_merged.rds")

cat("\nDone — CellChat objects saved\n")