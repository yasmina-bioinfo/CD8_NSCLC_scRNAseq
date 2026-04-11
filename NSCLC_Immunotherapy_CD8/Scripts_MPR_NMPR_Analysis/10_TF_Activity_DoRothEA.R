# =============================================================================
# Script 10 — TF Activity Inference with DoRothEA / decoupleR
# Dataset  : ImmunoT (GSE207422, Hu et al. 2023)
# Object   : 08_CD8_MPR_NMPR.rds
# Clusters : cell_state column
# Condition: MPR (n=10 patients, ~5,493 CD8 cells) vs NMPR (n=6, ~12,983 cells)
# Method   : mlm (multivariate linear model) — benchmark method per
#            Badia-i-Mompel et al. 2022, Bioinformatics Advances, vbac016
#            https://doi.org/10.1093/bioadv/vbac016
# Output   : Heatmap A (clusters × MPR/NMPR) + Heatmap B (clusters only)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------

library(Seurat)
library(decoupleR)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(tibble)
# -----------------------------------------------------------------------------
# 1. LOAD SEURAT OBJECT
# -----------------------------------------------------------------------------

seurat <- readRDS("objects/08_CD8_MPR_NMPR.rds")

# Verify key columns
stopifnot("cell_state" %in% colnames(seurat@meta.data))
stopifnot("PathResponse" %in% colnames(seurat@meta.data)) 

# Quick check
table(seurat$cell_state, seurat$PathResponse)

# -----------------------------------------------------------------------------
# 2. GET DoRothEA NETWORK FROM OmniPath
# -----------------------------------------------------------------------------
# DoRothEA confidence levels A-E; A = highest evidence
# We keep A-C for a balance of coverage and confidence
# Level D-E are highly inferred and noisier in single-cell context

net <- get_dorothea(organism = "human", levels = c("A", "B", "C"))

# Inspect network
head(net)
# Columns: source (TF), target (gene), mor (mode of regulation: +1/-1), weight

message("DoRothEA network: ", nrow(net), " interactions, ",
        length(unique(net$source)), " TFs")

# -----------------------------------------------------------------------------
# 3. PREPARE EXPRESSION MATRIX
# -----------------------------------------------------------------------------
# decoupleR works on the normalized log-counts matrix (RNA assay, data layer)

mat <- GetAssayData(seurat, assay = "RNA", layer = "data")
mat <- as.matrix(mat)

# -----------------------------------------------------------------------------
# 4. RUN TF ACTIVITY INFERENCE — METHOD: MLM
# -----------------------------------------------------------------------------
# mlm = multivariate linear model
# Models all targets of a TF simultaneously → controls inter-gene correlations
# Recommended benchmark method (Badia-i-Mompel et al. 2022)
# Output: activity score per cell per TF (z-score-like)

tf_acts <- run_mlm(
  mat       = mat,
  net       = net,
  .source   = "source",
  .target   = "target",
  .mor      = "mor",
  minsize   = 5          # minimum number of targets per TF (filters poorly covered TFs)
)

# tf_acts is a tibble: columns = statistic, source (TF), condition (cell), score, p_value
head(tf_acts)

# Keep only mlm statistic
tf_acts <- tf_acts %>% filter(statistic == "mlm")

message("TF activity computed for ", length(unique(tf_acts$source)), " TFs")

# -----------------------------------------------------------------------------
# 5. ADD TF ACTIVITY SCORES TO SEURAT METADATA
# -----------------------------------------------------------------------------
# Pivot to wide matrix: rows = cells, cols = TFs

tf_mat <- tf_acts %>%
  pivot_wider(id_cols = condition,
              names_from = source,
              values_from = score) %>%
  column_to_rownames("condition") %>%
  as.matrix()

# tf_mat: cells × TFs
# Add to Seurat as a new assay for convenient access

tf_assay <- CreateAssayObject(counts = t(tf_mat))
seurat[["TF"]] <- tf_assay
DefaultAssay(seurat) <- "TF"

# -----------------------------------------------------------------------------
# 6. SELECT TOP 25 TFs BY VARIANCE ACROSS CELLS
# -----------------------------------------------------------------------------

tf_var <- apply(tf_mat, 2, var, na.rm = TRUE)
top25_tfs <- names(sort(tf_var, decreasing = TRUE))[1:25]

message("Top 25 TFs by variance: ", paste(top25_tfs, collapse = ", "))

# -----------------------------------------------------------------------------
# 7. HEATMAP B — TF ACTIVITY BY CLUSTER (cell_state)
# Descriptive view: which TF dominates which CD8 state
# -----------------------------------------------------------------------------

# Compute mean TF activity per cluster
meta <- seurat@meta.data %>%
  select(cell_state) %>%
  rownames_to_column("cell")

tf_long <- tf_acts %>%
  filter(source %in% top25_tfs) %>%
  rename(cell = condition)

heatmap_B_data <- tf_long %>%
  left_join(meta, by = "cell") %>%
  group_by(source, cell_state) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = cell_state, values_from = mean_score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# Scale per TF (row) for visual contrast
heatmap_B_scaled <- t(scale(t(heatmap_B_data)))

# Color palette
col_palette <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
)(100)

# Plot & save
png("10_Heatmap_B_TF_by_Cluster.png", width = 2400, height = 2800,
    res = 450, bg = "white")

pheatmap(
  heatmap_B_scaled,
  color            = col_palette,
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  border_color     = NA,
  fontsize_row     = 7,
  fontsize_col     = 7,
  main             = "TF Activity by CD8 State\n(DoRothEA mlm, top 25 TFs by variance)",
  angle_col        = 45,
  na_col           = "grey90"
)

dev.off()
message("Heatmap B saved: 10_Heatmap_B_TF_by_Cluster.png")

# -----------------------------------------------------------------------------
# 8. HEATMAP A — TF ACTIVITY BY CLUSTER × CONDITION (MPR vs NMPR)
# Shows differential TF activity within each CD8 state between responders
# -----------------------------------------------------------------------------

meta2 <- seurat@meta.data %>%
  select(cell_state, PathResponse) %>%
  rownames_to_column("cell") %>%
  mutate(group = paste0(cell_state, "__", PathResponse))  # e.g. CD8_Exhausted_Terminal__MPR

heatmap_A_data <- tf_long %>%
  left_join(meta2, by = "cell") %>%
  group_by(source, group) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# Scale per TF (row)
heatmap_A_scaled <- t(scale(t(heatmap_A_data)))

# Column annotation: MPR vs NMPR
col_labels   <- colnames(heatmap_A_scaled)
condition_ann <- data.frame(
  Condition = ifelse(grepl("MPR$", col_labels) & !grepl("NMPR", col_labels),
                     "MPR", "NMPR"),
  row.names = col_labels
)
ann_colors <- list(
  Condition = c(MPR = "#2166AC", NMPR = "#D6604D")
)

# Clean column labels for display (remove redundant prefix)
display_names <- gsub("__MPR|__NMPR", "", col_labels)

png("10_Heatmap_A_TF_by_Cluster_Condition.png", width = 3600, height = 2800,
    res = 450, bg = "white")

pheatmap(
  heatmap_A_scaled,
  color             = col_palette,
  annotation_col    = condition_ann,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,   # Keep MPR/NMPR paired per cluster
  border_color      = NA,
  fontsize_row      = 7,
  fontsize_col      = 6,
  labels_col        = display_names,
  main              = "TF Activity by CD8 State × PathResponse (MPR vs NMPR)\n(DoRothEA mlm, top 25 TFs by variance)",
  angle_col         = 45,
  na_col            = "grey90",
  gaps_col          = seq(2, ncol(heatmap_A_scaled) - 1, by = 2)  # visual gap between clusters
)

dev.off()
message("Heatmap A saved: 10_Heatmap_A_TF_by_Cluster_Condition.png")

# -----------------------------------------------------------------------------
# 9. SESSION INFO (reproducibility)
# -----------------------------------------------------------------------------

sink("10_session_info.txt")
sessionInfo()
sink()

message("\n=== Script 10 complete ===")
message("Outputs:")
message("  10_Heatmap_B_TF_by_Cluster.png")
message("  10_Heatmap_A_TF_by_Cluster_Condition.png")
message("  10_UMAP_Top6_TF_Activity.png")
message("  10_session_info.txt")

