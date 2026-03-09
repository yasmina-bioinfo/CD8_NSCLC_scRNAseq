# ============================================================
# 01b_Add_PathResponse_TME.R
# Add Pathologic Response (MPR/NMPR) to TME Seurat object
# Source: GSE207422_NSCLC_scRNAseq_metadata.xlsx
# ============================================================

library(Seurat)

# Load TME object
seu_tme <- readRDS("Objects/04_TME_annotated.rds")

# Mapping from original GEO metadata
mapping <- data.frame(
  Sample = c("BD_immune01","BD_immune02","BD_immune03","BD_immune04","BD_immune05",
             "BD_immune06","BD_immune07","BD_immune08","BD_immune09","BD_immune10",
             "BD_immune11","BD_immune12","BD_immune13","BD_immune14","BD_immune15"),
  PathResponse = c("NE","NMPR","MPR","NMPR","NMPR",
                   "pCR","NMPR","NMPR","NMPR","NMPR",
                   "MPR","NMPR","NMPR","MPR","NMPR")
)

# Add to Seurat object
seu_tme$PathResponse <- mapping$PathResponse[match(seu_tme$Sample, mapping$Sample)]

# Verify
table(seu_tme$PathResponse)

# Save updated object
saveRDS(seu_tme, "Objects/04_TME_annotated.rds")
cat("Done — PathResponse added to TME\n")

# Subset to MPR and NMPR only (exclude pCR and NE)
seu_tme_mpr <- subset(seu_tme, PathResponse %in% c("MPR", "NMPR"))
table(seu_tme_mpr$PathResponse)

# Save
saveRDS(seu_tme_mpr, "Objects/04_TME_MPR_NMPR.rds")
cat("Done — TME MPR/NMPR object saved\n")