#!/usr/bin/env Rscript

# ============================================================
# 18_Merge_CD4_CD8_and_BarplotsSplit.R
# Merge annotated CD4 + CD8 objects and generate ONE figure:
# - Panel A: CD4_state proportions split by condition
# - Panel B: CD8_state proportions split by condition
# Uses a fixed, reusable palette (save it for Immunotherapy too).
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# Inputs
# -----------------------------
cd4_rds <- "Objects/16_CD4_clean_annotated.rds"
cd8_rds <- "Objects/17_CD8_clean_annotated.rds"

# -----------------------------
# Output
# -----------------------------
fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(fig_dir, "18_CD4_CD8_states_barplot_splitByCondition.png")
out_pdf <- file.path(fig_dir, "18_CD4_CD8_states_barplot_splitByCondition.pdf")

# ============================================================
# Reproducible palette (REUSE this exact block in Immunotherapy)
# ============================================================
state_cols <- c(
  "CD4_Naive_CM"          = "#2166AC",
  "CD4_Effector_Th1_like" = "#F97B2B",
  "Treg_Activated"        = "#F4A7C3",
  
  "CD8_Naive_CM"       = "#A8D8E8",  
  "CD8_Effector_GZMK"  = "#FDBB84",  
  "CD8_TRM_Cytotoxic"  = "#D73027",  
  "CD8_Proliferating"  = "#CBC9E2"   
)

# Optional: define stable order (useful for legend + stacks)
state_levels <- c(
  "CD4_Naive_CM", "CD4_Effector_Th1_like", "Treg_Activated",
  "CD8_Naive_CM", "CD8_Effector_GZMK", "CD8_TRM_Cytotoxic", "CD8_Proliferating"
)

# -----------------------------
# Load
# -----------------------------
seu_cd4 <- readRDS(cd4_rds)
seu_cd8 <- readRDS(cd8_rds)

# Guardrails
stopifnot("Sample_Origin" %in% colnames(seu_cd4@meta.data))
stopifnot("Sample_Origin" %in% colnames(seu_cd8@meta.data))
stopifnot("CD4_state" %in% colnames(seu_cd4@meta.data))
stopifnot("CD8_state" %in% colnames(seu_cd8@meta.data))

# Drop empty levels (clean tables)
seu_cd4$Sample_Origin <- droplevels(seu_cd4$Sample_Origin)
seu_cd8$Sample_Origin <- droplevels(seu_cd8$Sample_Origin)

# -----------------------------
# Merge
# -----------------------------
# Make sure cell names are unique across objects
seu_cd4 <- RenameCells(seu_cd4, add.cell.id = "CD4")
seu_cd8 <- RenameCells(seu_cd8, add.cell.id = "CD8")

seu_T <- merge(seu_cd4, y = seu_cd8, project = "LUAD_Tcells_CD4_CD8")

# -----------------------------
# Build a tidy table for plotting (NO additional figure before this)
# -----------------------------
meta <- seu_T[[]] %>%
  transmute(
    Sample_Origin,
    CD4_state = ifelse(is.na(CD4_state), NA, as.character(CD4_state)),
    CD8_state = ifelse(is.na(CD8_state), NA, as.character(CD8_state))
  )

# Long format: one column "panel" (CD4 vs CD8), one "state"
df_cd4 <- meta %>%
  filter(!is.na(CD4_state)) %>%
  transmute(panel = "CD4", Sample_Origin, state = CD4_state)

df_cd8 <- meta %>%
  filter(!is.na(CD8_state)) %>%
  transmute(panel = "CD8", Sample_Origin, state = CD8_state)

df <- bind_rows(df_cd4, df_cd8) %>%
  mutate(
    state = factor(state, levels = state_levels),
    panel = factor(panel, levels = c("CD4", "CD8"))
  ) %>%
  count(panel, Sample_Origin, state, name = "n") %>%
  group_by(panel, Sample_Origin) %>%
  mutate(Percent = 100 * n / sum(n)) %>%
  ungroup()

# -----------------------------
# One figure: facets (CD4 + CD8) split by condition
# -----------------------------
p <- ggplot(df, aes(x = Sample_Origin, y = Percent, fill = state)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  theme_bw() +
  ylab("Percentage within panel (CD4 or CD8)") +
  xlab("") +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text  = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(5, "mm"),
    legend.spacing.y = unit(0.4, "cm")
  )

cd8_only_cols <- state_cols[names(state_cols) %in% 
                              c("CD8_Naive_CM", "CD8_Effector_GZMK", 
                                "CD8_TRM_Cytotoxic", "CD8_Proliferating")]

p_cd8_only <- ggplot(df %>% filter(panel == "CD8"),
                     aes(x = Sample_Origin, y = Percent, fill = state)) +
  geom_bar(stat = "identity", width = 0.85) +
  scale_fill_manual(values = cd8_only_cols, drop = TRUE) +
  theme_bw() +
  ylab("Percentage of cells") +
  xlab("") +
  theme(
    axis.text.x  = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 12),
    legend.title = element_blank(),
    legend.text  = element_text(size = 13),
    legend.key.size = unit(6, "mm"),
    legend.spacing.y = unit(0.6, "cm")
  )

ggsave(file.path(fig_dir, "18_CD8_states_barplot_splitByCondition.png"),
       p_cd8_only, width = 7, height = 6, dpi = 450, bg = "white")


ggsave(out_png, plot = p, width = 8, height = 8, dpi = 450)
ggsave(out_pdf, plot = p, width = 8, height = 8)

# -----------------------------
# Save merged object (optional but useful)
# -----------------------------
saveRDS(seu_T, "Objects/18_Tcells_CD4_CD8_merged.rds")

cat("\nSaved figure -> ", out_png, "\n", sep = "")
cat("Saved merged object -> Objects/18_Tcells_CD4_CD8_merged.rds\n")
cat("\nPalette to reuse (state_cols):\n")
print(state_cols)
