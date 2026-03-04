# ============================================================
# 23b_ProportionTest_ImmunoT.R
# Chi-square test sur les proportions CD8 — PR vs SD
# Tous les clusters CD8
# ============================================================

library(ggplot2)
library(tidyr)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Contingency table — tous les clusters
# -----------------------------
counts <- matrix(
  c(1236,   16,   # CD8_Activated_HLAII_high
    1054,  755,   # CD8_Early_Activated_NR4A_high
    4914,  767,   # CD8_Effector_GZMK
    3741,   41,   # CD8_Exhausted_Terminal
    27,  860,   # CD8_IFN_Stress_Response
    1866,  375,   # CD8_Proliferating
    1465, 1007,   # CD8_Terminal_CX3CR1
    1667,  528),  # CD8_TRM_like
  nrow = 8, byrow = TRUE,
  dimnames = list(
    c("CD8_Activated_HLAII_high", "CD8_Early_Activated_NR4A_high",
      "CD8_Effector_GZMK", "CD8_Exhausted_Terminal",
      "CD8_IFN_Stress_Response", "CD8_Proliferating",
      "CD8_Terminal_CX3CR1", "CD8_TRM_like"),
    c("PR", "SD")
  )
)

# -----------------------------
# Global chi-square test
# -----------------------------
chi_test <- chisq.test(counts)
cat("Chi-square test:\n")
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
    PR      = counts[cluster, "PR"],
    SD      = counts[cluster, "SD"],
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
# Colors — palette ImmunoT
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
df_plot$cluster <- rownames(counts)
df_plot$prop_PR <- df_plot$PR / total_n["PR"] * 100
df_plot$prop_SD <- df_plot$SD / total_n["SD"] * 100

df_long <- pivot_longer(df_plot,
                        cols      = c("prop_PR", "prop_SD"),
                        names_to  = "condition",
                        values_to = "proportion")
df_long$condition <- ifelse(df_long$condition == "prop_PR", "PR", "SD")
df_long$condition <- factor(df_long$condition, levels = c("PR", "SD"))
df_long <- merge(df_long, results[, c("cluster", "sig")], by = "cluster")

p <- ggplot(df_long, aes(x = condition, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 12, color = "#000000"),
    axis.text.y      = element_text(size = 11, color = "#000000"),
    legend.title     = element_text(size = 11),
    plot.title       = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle    = element_text(size = 11, color = "#555555", face = "italic"),
    plot.caption     = element_text(size = 11, color = "#555555")
  ) +
  labs(
    x        = "",
    y        = "Percentage of CD8 T cells (%)",
    fill     = "CD8 state",
    title    = "CD8 T cell state proportions — PR vs SD",
    subtitle = "Fisher's exact test with BH correction",
    caption  = "Fisher's exact test, BH correction. ***p<0.001"
  )

# Star positions (proportion) = visual midpoint of each colored zone on the barplot
# To adjust: look at the barplot, estimate the center of the cluster's zone
# Example: if CD8_Exhausted_Terminal spans 38%-57% in PR, set proportion = 47
# condition = the condition where the cluster is enriched (OR > 1 = PR, OR < 1 = SD)
sig_data <- data.frame(
  cluster    = c("CD8_Activated_HLAII_high", "CD8_Early_Activated_NR4A_high",
                 "CD8_Effector_GZMK", "CD8_Exhausted_Terminal",
                 "CD8_Proliferating", "CD8_TRM_like",
                 "CD8_Early_Activated_NR4A_high", "CD8_IFN_Stress_Response",
                 "CD8_Terminal_CX3CR1"),
  condition  = c("PR", "PR", "PR", "PR", "PR", "PR", "SD", "SD", "SD"),
  proportion = c(3.9, 13, 29.7, 56.8, 74.6, 85.6, 16, 46.2, 88.4),
  sig        = c("***", "***", "***", "***", "***", "**", "***", "***", "***")
)

p <- p + geom_text(
  data        = sig_data,
  aes(x = condition, y = proportion, label = sig),
  inherit.aes = FALSE,
  size        = 5,
  color       = "#000000",
  vjust       = -0.5
)


ggsave(file.path(fig_dir, "19_ProportionTest_ImmunoT.png"),
       p, width = 8, height = 6, dpi = 450, bg = "white")

cat("\nDone\n")