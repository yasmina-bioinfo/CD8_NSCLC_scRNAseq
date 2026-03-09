# ============================================================
# 05_ProportionTest_CD8_MPR_NMPR.R
# Fisher's exact test sur les proportions CD8 — MPR vs NMPR
# ============================================================

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

seu_mpr <- readRDS("Objects/08_CD8_MPR_NMPR.rds")

# -----------------------------
# Contingency table — calculée dynamiquement
# -----------------------------
counts_df <- seu_mpr@meta.data %>%
  count(cell_state, PathResponse) %>%
  pivot_wider(names_from = PathResponse, values_from = n, values_fill = 0)

counts <- as.matrix(counts_df[, c("MPR", "NMPR")])
rownames(counts) <- counts_df$cell_state

cat("Contingency table:\n")
print(counts)

# -----------------------------
# Global chi-square test
# -----------------------------
chi_test <- chisq.test(counts)
cat("\nChi-square test:\n")
print(chi_test)

# -----------------------------
# Pairwise Fisher tests per cluster
# -----------------------------
results <- data.frame()

for (cluster in rownames(counts)) {
  in_cluster  <- counts[cluster, ]
  out_cluster <- colSums(counts) - in_cluster
  mat <- rbind(in_cluster, out_cluster)
  ft  <- fisher.test(mat)
  results <- rbind(results, data.frame(
    cluster = cluster,
    MPR     = counts[cluster, "MPR"],
    NMPR    = counts[cluster, "NMPR"],
    OR      = round(ft$estimate, 3),
    p_value = round(ft$p.value, 4),
    p_adj   = NA
  ))
}

results$p_adj <- round(p.adjust(results$p_value, method = "BH"), 4)
results$sig   <- ifelse(results$p_adj < 0.001, "***",
                        ifelse(results$p_adj < 0.01,  "**",
                               ifelse(results$p_adj < 0.05,  "*", "ns")))

cat("\nPairwise Fisher tests:\n")
print(results)

# -----------------------------
# Colors
# -----------------------------
cluster_colors <- c(
  "CD8_Effector_GZMK"             = "#92C5DE",
  "CD8_Exhausted_Terminal"        = "#D73027",
  "CD8_Terminal_CX3CR1"           = "#FDBB84",
  "CD8_Proliferating"             = "#CBC9E2",
  "CD8_TRM_like"                  = "#D6A5A5",
  "CD8_Early_Activated_NR4A_high" = "#A1D99B",
  "CD8_Activated_HLAII_high"      = "#D9D97A",
  "CD8_IFN_Stress_Response"       = "#A8D8E8"
)

# -----------------------------
# Plot
# -----------------------------
total_n <- colSums(counts)
df_plot <- as.data.frame(counts)
df_plot$cluster  <- rownames(counts)
df_plot$prop_MPR  <- df_plot$MPR  / total_n["MPR"]  * 100
df_plot$prop_NMPR <- df_plot$NMPR / total_n["NMPR"] * 100

df_long <- pivot_longer(df_plot,
                        cols      = c("prop_MPR", "prop_NMPR"),
                        names_to  = "condition",
                        values_to = "proportion")
df_long$condition <- ifelse(df_long$condition == "prop_MPR", "MPR", "NMPR")
df_long$condition <- factor(df_long$condition, levels = c("MPR", "NMPR"))
df_long <- merge(df_long, results[, c("cluster", "sig")], by = "cluster")

p <- ggplot(df_long, aes(x = condition, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x   = element_text(size = 12, color = "#000000"),
    axis.text.y   = element_text(size = 11, color = "#000000"),
    legend.title  = element_text(size = 11),
    legend.text   = element_text(size = 11),
    plot.title    = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle = element_text(size = 11, color = "#555555", face = "italic"),
    plot.caption  = element_text(size = 11, color = "#555555")
  ) +
  labs(
    x        = "",
    y        = "Percentage of CD8 T cells (%)",
    fill     = "CD8 state",
    title    = "CD8 T cell state proportions — MPR vs NMPR",
    subtitle = "Fisher's exact test with BH correction",
    caption  = "Fisher's exact test, BH correction. ***p<0.001"
  )

sig_data <- data.frame(
  cluster   = c("CD8_TRM_like", "CD8_Terminal_CX3CR1", "CD8_Proliferating",
                "CD8_Exhausted_Terminal", "CD8_Effector_GZMK",
                "CD8_IFN_Stress_Response", "CD8_Early_Activated_NR4A_high",
                "CD8_Activated_HLAII_high"),
  condition = c("MPR", "MPR", "MPR", "MPR", "MPR",
                "NMPR", "NMPR", "NMPR"),
  proportion = c(3, 12, 23, 45, 82, 42, 84, 95),
  sig        = c("*", "*", "***", "***", "***", "***", "***", "***")
)

p <- p + geom_text(
  data        = sig_data,
  aes(x = condition, y = proportion, label = sig),
  inherit.aes = FALSE,
  size        = 5,
  color       = "#000000",
  vjust       = 0.5
)

ggsave(file.path(fig_dir, "05_ProportionTest_CD8_MPR_NMPR.png"),
       p, width = 10, height = 6, dpi = 450, bg = "white")

cat("\nDone\n")