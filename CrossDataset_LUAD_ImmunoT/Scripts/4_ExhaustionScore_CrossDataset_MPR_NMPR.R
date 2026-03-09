#!/usr/bin/env Rscript
# ============================================================
# 08_ExhaustionScore_CrossDataset_MPR_NMPR.R
# Exhaustion Terminal Module Score — LUAD vs ImmunoT
# Cross-dataset comparison — nLung/tLung + MPR/NMPR
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
seu_luad   <- readRDS("Objects/15_CD8_annotated_states.rds")
seu_immuno <- readRDS("Objects/08_CD8_MPR_NMPR.rds")

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
  condition = seu_immuno$PathResponse,
  dataset   = "ImmunoT (GSE207422)"
)

df_all <- bind_rows(df_luad, df_immuno)
df_all$condition <- factor(df_all$condition,
                           levels = c("nLung", "tLung", "MPR", "NMPR"))

# -----------------------------
# Wilcoxon tests
# -----------------------------
wt_luad <- wilcox.test(
  df_luad$score[df_luad$condition == "nLung"],
  df_luad$score[df_luad$condition == "tLung"]
)
wt_immuno <- wilcox.test(
  df_immuno$score[df_immuno$condition == "MPR"],
  df_immuno$score[df_immuno$condition == "NMPR"]
)

format_pval <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(paste0("p = ", round(p, 3)))
  return(paste0("p = ", round(p, 2)))
}

label_luad   <- format_pval(wt_luad$p.value)
label_immuno <- format_pval(wt_immuno$p.value)

cat("Wilcoxon LUAD (nLung vs tLung):", label_luad, "\n")
cat("Wilcoxon ImmunoT (MPR vs NMPR):", label_immuno, "\n")

# -----------------------------
# Annotation dataframe for p-values
# -----------------------------
annot_df <- data.frame(
  dataset   = c("LUAD (GSE131907)", "ImmunoT (GSE207422)"),
  label     = c(label_luad, label_immuno),
  x         = c(1.5, 1.5),
  y         = c(
    max(df_luad$score)   * 1.05,
    max(df_immuno$score) * 1.05
  )
)

# -----------------------------
# Colors
# -----------------------------
condition_colors <- c(
  "nLung" = "#4575B4",
  "tLung" = "#D73027",
  "MPR"   = "#4575B4",
  "NMPR"  = "#D73027"
)

# -----------------------------
# Plot — Violin + Boxplot + Wilcoxon
# -----------------------------
p <- ggplot(df_all, aes(x = condition, y = score, fill = condition)) +
  geom_violin(alpha = 0.8, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  geom_text(
    data    = annot_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4, fontface = "italic", color = "#333333"
  ) +
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

ggsave(file.path(fig_dir, "08_ExhaustionScore_CrossDataset_MPR_NMPR.png"),
       p, width = 10, height = 6, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "08_ExhaustionScore_CrossDataset_MPR_NMPR.pdf"),
       p, width = 10, height = 6, bg = "white")

cat("\nDone — ExhaustionScore CrossDataset MPR/NMPR saved\n")
