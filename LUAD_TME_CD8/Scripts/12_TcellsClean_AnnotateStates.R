#!/usr/bin/env Rscript

# ============================================================
# 12_TcellsClean_AnnotateStates.R
# Annotate reclustered clean T cells with biological states
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# Input / Output
# -----------------------------
in_rds  <- "Objects/10_Tcells_clean_reclustered.rds"
out_rds <- "Objects/12_Tcells_clean_annotated.rds"

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load
# -----------------------------
seu <- readRDS(in_rds)

# Guardrails (fail fast if wrong object)
stopifnot("seurat_clusters" %in% colnames(seu@meta.data))
stopifnot("Sample_Origin" %in% colnames(seu@meta.data))

Idents(seu) <- "seurat_clusters"

# -----------------------------
# Cluster → State mapping (ROBUST)
# -----------------------------
cluster_map <- c(
  "0" = "T_Naive_CM",
  "1" = "NK_like_Cytotoxic",
  "2" = "T_TRM_like",
  "3" = "CD8_Effector_GZMK",
  "4" = "CD8_TRM",
  "5" = "T_Activated_Ig+",
  "6" = "Treg_Activated",
  "7" = "CD8_TRM_CXCL13",
  "8" = "T_Proliferating"
)

cluster_id <- as.character(seu$seurat_clusters)
seu$Tcell_state <- unname(cluster_map[cluster_id])

stopifnot("Tcell_state" %in% colnames(seu@meta.data))
stopifnot(!any(is.na(seu$Tcell_state)))

# -----------------------------
# Order + Color palette (define AFTER Tcell_state exists)
# -----------------------------
state_levels <- c(
  "T_Naive_CM",
  "T_TRM_like",
  "CD8_TRM",
  "CD8_TRM_CXCL13",
  "CD8_Effector_GZMK",
  "Treg_Activated",
  "T_Activated_Ig+",
  "T_Proliferating",
  "NK_like_Cytotoxic"
)

seu$Tcell_state <- factor(seu$Tcell_state, levels = state_levels)
Idents(seu) <- "Tcell_state"

state_cols <- c(
  "T_Naive_CM"        = "#1F77B4",
  "T_TRM_like"        = "#9467BD",
  "CD8_TRM"           = "#FF7F0E",
  "CD8_Effector_GZMK" = "#D62728",
  "CD8_TRM_CXCL13"    = "#2CA02C",
  "Treg_Activated"    = "#E377C2",
  "T_Proliferating"   = "#17BECF",
  "T_Activated_Ig+"   = "#8C564B",
  "NK_like_Cytotoxic" = "#7F7F7F"
)

# -----------------------------
# UMAP global (by state)
# -----------------------------
p1 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "Tcell_state",
  cols = state_cols,
  label = TRUE
) + theme_bw()

ggsave(file.path(fig_dir, "12_TcellsClean_UMAP_states.png"),
       plot = p1, width = 7, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "12_TcellsClean_UMAP_states.pdf"),
       plot = p1, width = 7, height = 6)

# -----------------------------
# UMAP split by condition
# -----------------------------
p2 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "Tcell_state",
  split.by = "Sample_Origin",
  cols = state_cols,
  label = TRUE,
  ncol = 2
) + theme_bw()

ggsave(file.path(fig_dir, "12_TcellsClean_UMAP_states_split.png"),
       plot = p2, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "12_TcellsClean_UMAP_states_split.pdf"),
       plot = p2, width = 12, height = 6)

# -----------------------------
# Barplot – Proportion of T cell states by condition
# -----------------------------
df <- seu[[]] %>%
  dplyr::count(Sample_Origin, Tcell_state, name = "n") %>%
  dplyr::group_by(Sample_Origin) %>%
  dplyr::mutate(Percent = 100 * n / sum(n)) %>%
  dplyr::ungroup()

p_bar <- ggplot(df, aes(x = Sample_Origin, y = Percent, fill = Tcell_state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  theme_bw() +
  ylab("Percentage of T cells") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "12_TcellsClean_barplot_states_by_condition.png"),
       plot = p_bar, width = 7, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "12_TcellsClean_barplot_states_by_condition.pdf"),
       plot = p_bar, width = 7, height = 6)

# -----------------------------
# DotPlot – validation markers per state (clean)
# -----------------------------

markers <- c(
  # Naive / CM
  "IL7R", "CCR7", "TCF7", "SELL",
  
  # CD8 cytotoxic / effector
  "CD8A", "CD8B", "GZMK", "NKG7", "PRF1", "GZMB",
  
  # TRM
  "ZNF683", "ITGAE", "CXCR6",
  
  # Treg
  "FOXP3", "IL2RA", "CTLA4", "TIGIT",
  
  # CXCL13 TIL
  "CXCL13",
  
  # Proliferation
  "MKI67", "TOP2A"
)

p_dot <- DotPlot(seu, features = markers) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    )
  )

ggsave(file.path(fig_dir, "12_TcellsClean_DotPlot_states.png"),
       plot = p_dot, width = 12, height = 6, dpi = 300)

ggsave(file.path(fig_dir, "12_TcellsClean_DotPlot_states.pdf"),
       plot = p_dot, width = 12, height = 6)

# -----------------------------
# Validation + Save
# -----------------------------
cat("\n=== 12_TcellsClean_AnnotateStates : VALIDATION ===\n")
cat("Cells per state:\n")
print(table(seu$Tcell_state))
cat("\nCells per condition:\n")
print(table(seu$Sample_Origin))

saveRDS(seu, out_rds)
cat("\nSaved -> ", out_rds, "\n")