# ============================================================
# 01_Add_PathResponse_metadata.R
# Add Pathologic Response (MPR/NMPR) to CD8 Seurat object
# Source: GSE207422_NSCLC_scRNAseq_metadata.xlsx
# ============================================================

library(Seurat)

# Load CD8 object
seu_cd8 <- readRDS("Objects/07_CD8_annotated.rds")

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
seu_cd8$PathResponse <- mapping$PathResponse[match(seu_cd8$Sample, mapping$Sample)]

# Verify
table(seu_cd8$PathResponse)

# Save updated object
saveRDS(seu_cd8, "Objects/07_CD8_annotated.rds")

cat("Done — PathResponse added\n")

# Subset to MPR and NMPR only
seu_mpr <- subset(seu_cd8, PathResponse %in% c("MPR", "NMPR"))
table(seu_mpr$PathResponse)

# Save
saveRDS(seu_mpr, "Objects/08_CD8_MPR_NMPR.rds")
cat("Done — MPR/NMPR object saved\n")