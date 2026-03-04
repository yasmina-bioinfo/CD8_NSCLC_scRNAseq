# ============================================================
# 21_ProportionTest_LUAD.R
# Chi-square test sur les proportions CD8 — nLung vs tLung
# ============================================================

library(ggplot2)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Contingency table
# -----------------------------
counts <- matrix(
  c(132, 252,   # CD8_Effector_GZMK
    245, 243,   # CD8_Naive_CM
    13,  34,   # CD8_Proliferating
    17, 223),  # CD8_TRM_Cytotoxic
  nrow = 4, byrow = TRUE,
  dimnames = list(
    c("CD8_Effector_GZMK", "CD8_Naive_CM", 
      "CD8_Proliferating", "CD8_TRM_Cytotoxic"),
    c("nLung", "tLung")
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
  # 2x2 table: cluster vs rest
  in_cluster  <- counts[cluster, ]
  out_cluster <- colSums(counts) - in_cluster
  
  mat <- rbind(in_cluster, out_cluster)
  ft  <- fisher.test(mat)
  
  results <- rbind(results, data.frame(
    cluster  = cluster,
    nLung    = counts[cluster, "nLung"],
    tLung    = counts[cluster, "tLung"],
    OR       = round(ft$estimate, 3),
    p_value  = round(ft$p.value, 4),
    p_adj    = NA
  ))
}

# FDR correction
results$p_adj <- round(p.adjust(results$p_value, method = "BH"), 4)
results$sig   <- ifelse(results$p_adj < 0.001, "***",
                        ifelse(results$p_adj < 0.01,  "**",
                               ifelse(results$p_adj < 0.05,  "*", "ns")))

cat("\nPairwise Fisher tests:\n")
print(results)

# -----------------------------
# Plot — proportions + significance
# -----------------------------
total_n <- colSums(counts)
df_plot <- as.data.frame(counts)
df_plot$cluster   <- rownames(counts)
df_plot$prop_nLung <- df_plot$nLung / total_n["nLung"] * 100
df_plot$prop_tLung <- df_plot$tLung / total_n["tLung"] * 100

df_long <- tidyr::pivot_longer(df_plot, 
                               cols = c("prop_nLung", "prop_tLung"),
                               names_to = "condition",
                               values_to = "proportion")
df_long$condition <- ifelse(df_long$condition == "prop_nLung", "nLung", "tLung")
df_long$condition <- factor(df_long$condition, levels = c("nLung", "tLung"))

# Add significance
df_long <- merge(df_long, results[, c("cluster", "sig")], by = "cluster")

cluster_colors <- c(
  "CD8_Naive_CM"       = "#A8D8E8",
  "CD8_Effector_GZMK"  = "#FDBB84",
  "CD8_TRM_Cytotoxic"  = "#D73027",
  "CD8_Proliferating"  = "#CBC9E2"
)

p <- ggplot(df_long, aes(x = condition, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 12, color = "#000000"),
    axis.text.y      = element_text(size = 11, color = "#000000"),
    legend.title     = element_text(size = 11),
    plot.title       = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle    = element_text(size = 11, color = "#555555", face = "italic")
  ) +
  labs(
    x        = "",
    y        = "Percentage of CD8 T cells (%)",
    fill     = "CD8 state",
    title    = "CD8 T cell state proportions — nLung vs tLung",
    subtitle = "Fisher's exact test with BH correction",
    caption = "Fisher's exact test, BH correction. ***p<0.001"
  )

# Significance annotations
sig_data <- data.frame(
  cluster    = c("CD8_Naive_CM", "CD8_TRM_Cytotoxic"),
  condition  = c("nLung", "tLung"),
  proportion = c(50, 15), # this is not a biological value mais "*" position
  sig        = c("***", "***")
)

p <- p + geom_text(
  data        = sig_data,
  aes(x = condition, y = proportion, label = sig),
  inherit.aes = FALSE,
  size        = 5,
  color       = "#000000",
  vjust       = -0.5
)

ggsave(file.path(fig_dir, "23a_ProportionTest_LUAD.png"),
       p, width = 8, height = 6, dpi = 450, bg = "white")

cat("\nDone\n")