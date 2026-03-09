#!/usr/bin/env Rscript
# ============================================================
# 2_ExhaustionScore_CrossDataset.R
# Exhaustion Terminal Module Score — LUAD vs ImmunoT
# Cross-dataset comparison
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Exhaustion signature
# -----------------------------
exhaustion_genes <- c("PDCD1", "TIGIT", "HAVCR2", "LAG3", 
                      "CTLA4", "ENTPD1", "LAYN", "PRDM1")

# -----------------------------
# Load objects
# -----------------------------

seu_luad   <- readRDS("../LUAD/Objects/15_CD8_annotated_states.rds")
seu_immuno <- readRDS("../Immunotherapy/Objects/07_CD8_annotated.rds")

# -----------------------------
# Calculate module score
# -----------------------------
seu_luad <- AddModuleScore(
  seu_luad,
  features = list(exhaustion_genes),
  name     = "Exhaustion_Terminal"
)

seu_immuno <- AddModuleScore(
  seu_immuno,
  features = list(exhaustion_genes),
  name     = "Exhaustion_Terminal"
)

# -----------------------------
# Extract scores
# -----------------------------
df_luad <- data.frame(
  score     = seu_luad$Exhaustion_Terminal1,
  condition = seu_luad$Sample_Origin,
  dataset   = "LUAD (GSE131907)"
)

df_immuno <- data.frame(
  score     = seu_immuno$Exhaustion_Terminal1,
  condition = seu_immuno$RECIST,
  dataset   = "ImmunoT (GSE207422)"
)

df_all <- bind_rows(df_luad, df_immuno)
df_all$condition <- factor(df_all$condition, 
                           levels = c("nLung", "tLung", "PR", "SD"))

# -----------------------------
# Colors
# -----------------------------
condition_colors <- c(
  "nLung" = "#4575B4",
  "tLung" = "#D73027",
  "PR"    = "#4575B4",
  "SD"    = "#D73027"
)

# -----------------------------
# Plot — Violin + Boxplot
# -----------------------------
p <- ggplot(df_all, aes(x = condition, y = score, fill = condition)) +
  geom_violin(alpha = 0.8, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = condition_colors) +
  facet_wrap(~ dataset, scales = "free_x") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 12, color = "#000000"),
    axis.text.y      = element_text(size = 11, color = "#000000"),
    strip.text       = element_text(size = 12, face = "bold", color = "#000000"),
    strip.background = element_rect(fill = "#F5F5F5"),
    legend.position  = "none",
    plot.title       = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle    = element_text(size = 11, color = "#555555", face = "italic"),
    plot.caption     = element_text(size = 11, color = "#555555", face = "italic")
  ) +
  labs(
    x        = "",
    y        = "Exhaustion Terminal Module Score",
    title    = "Exhaustion Terminal Score — cross-dataset comparison",
    subtitle = "PDCD1, TIGIT, HAVCR2, LAG3, CTLA4, ENTPD1, LAYN, PRDM1",
    caption  = "Note: datasets are independent; comparison is qualitative only."
  )

ggsave(file.path(fig_dir, "ExhaustionScore_CrossDataset.png"),
       p, width = 10, height = 6, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "ExhaustionScore_CrossDataset.pdf"),
       p, width = 10, height = 6, bg = "white")

cat("\nDone\n")