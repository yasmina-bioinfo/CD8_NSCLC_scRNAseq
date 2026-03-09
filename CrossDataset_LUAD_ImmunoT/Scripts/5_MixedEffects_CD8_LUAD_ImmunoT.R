#!/usr/bin/env Rscript
# ============================================================
# 5_MixedEffects_CD8_LUAD_ImmunoT.R
# Mixed-effects modeling — CD8 cluster proportions + ModuleScore
# LUAD (nLung vs tLung) + ImmunoT (MPR vs NMPR)
# Patient as random effect to account for pseudoreplication
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

fig_dir    <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/CrossDataset/Results/figures"
output_dir <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/CrossDataset/Results/tables"
dir.create(fig_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Exhaustion signature
# -----------------------------
exhaustion_genes <- c("PDCD1", "TIGIT", "HAVCR2", "LAG3",
                      "CTLA4", "ENTPD1", "LAYN", "PRDM1")

# ============================================================
# PART 1 — LUAD (nLung vs tLung)
# ============================================================
seu_luad <- readRDS("C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/LUAD/Objects/16_CD8_LUAD_with_Patient.rds")
DefaultAssay(seu_luad) <- "RNA"

# Filter to nLung and tLung only
seu_luad <- subset(seu_luad, Sample_Origin %in% c("nLung", "tLung"))

# Add ModuleScore
seu_luad <- AddModuleScore(
  seu_luad,
  features = list(exhaustion_genes),
  name     = "Exhaustion_Terminal"
)

# -----------------------------
# Proportions per patient per cluster — LUAD
# -----------------------------
meta_luad <- seu_luad@meta.data %>%
  group_by(Patient, Sample_Origin, CD8_state) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(Patient, Sample_Origin) %>%
  mutate(total = sum(n_cells),
         proportion = n_cells / total) %>%
  ungroup() %>%
  mutate(Sample_Origin = factor(Sample_Origin, levels = c("nLung", "tLung")))

# -----------------------------
# Mixed-effects model per cluster — LUAD proportions
# -----------------------------
clusters_luad <- unique(meta_luad$CD8_state)
results_luad_prop <- lapply(clusters_luad, function(cl) {
  df <- meta_luad %>% filter(CD8_state == cl)
  # Need at least 2 patients per condition
  if (nrow(df) < 4) return(NULL)
  tryCatch({
    m <- lmer(proportion ~ Sample_Origin + (1 | Patient), data = df)
    coef_df <- as.data.frame(summary(m)$coefficients)
    data.frame(
      cluster   = cl,
      estimate  = coef_df["Sample_OriginTLung", "Estimate"],
      se        = coef_df["Sample_OriginTLung", "Std. Error"],
      p_value   = coef_df["Sample_OriginTLung", "Pr(>|t|)"],
      dataset   = "LUAD"
    )
  }, error = function(e) NULL)
})
results_luad_prop <- bind_rows(results_luad_prop)
results_luad_prop$p_adj <- p.adjust(results_luad_prop$p_value, method = "BH")

# -----------------------------
# Mixed-effects model — LUAD ModuleScore
# -----------------------------
meta_luad_score <- seu_luad@meta.data %>%
  group_by(Patient, Sample_Origin) %>%
  summarise(mean_score = mean(Exhaustion_Terminal1), .groups = "drop") %>%
  mutate(Sample_Origin = factor(Sample_Origin, levels = c("nLung", "tLung")))

m_luad_score <- lmer(mean_score ~ Sample_Origin + (1 | Patient),
                     data = meta_luad_score)
luad_score_result <- as.data.frame(summary(m_luad_score)$coefficients)

# ============================================================
# PART 2 — ImmunoT (MPR vs NMPR)
# ============================================================
seu_immuno <- readRDS("C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/Objects/08_CD8_MPR_NMPR.rds")
DefaultAssay(seu_immuno) <- "RNA"

# Add ModuleScore
seu_immuno <- AddModuleScore(
  seu_immuno,
  features = list(exhaustion_genes),
  name     = "Exhaustion_Terminal"
)

# -----------------------------
# Proportions per patient per cluster — ImmunoT
# -----------------------------
meta_immuno <- seu_immuno@meta.data %>%
  rename(Patient = Sample) %>%
  group_by(Patient, PathResponse, cell_state) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(Patient, PathResponse) %>%
  mutate(total = sum(n_cells),
         proportion = n_cells / total) %>%
  ungroup() %>%
  mutate(PathResponse = factor(PathResponse, levels = c("NMPR", "MPR")))

# -----------------------------
# Mixed-effects model per cluster — ImmunoT proportions
# -----------------------------
clusters_immuno <- unique(meta_immuno$cell_state)
results_immuno_prop <- lapply(clusters_immuno, function(cl) {
  df <- meta_immuno %>% filter(cell_state == cl)
  if (nrow(df) < 4) return(NULL)
  tryCatch({
    m <- lmer(proportion ~ PathResponse + (1 | Patient), data = df)
    coef_df <- as.data.frame(summary(m)$coefficients)
    data.frame(
      cluster  = cl,
      estimate = coef_df["PathResponseMPR", "Estimate"],
      se       = coef_df["PathResponseMPR", "Std. Error"],
      p_value  = coef_df["PathResponseMPR", "Pr(>|t|)"],
      dataset  = "ImmunoT"
    )
  }, error = function(e) NULL)
})
results_immuno_prop <- bind_rows(results_immuno_prop)
results_immuno_prop$p_adj <- p.adjust(results_immuno_prop$p_value, method = "BH")

# -----------------------------
# Mixed-effects model — ImmunoT ModuleScore
# -----------------------------
meta_immuno_score <- seu_immuno@meta.data %>%
  rename(Patient = Sample) %>%
  group_by(Patient, PathResponse) %>%
  summarise(mean_score = mean(Exhaustion_Terminal1), .groups = "drop") %>%
  mutate(PathResponse = factor(PathResponse, levels = c("NMPR", "MPR")))

m_immuno_score <- lmer(mean_score ~ PathResponse + (1 | Patient),
                       data = meta_immuno_score)
immuno_score_result <- as.data.frame(summary(m_immuno_score)$coefficients)

# ============================================================
# RESULTS SUMMARY
# ============================================================
cat("\n=== LUAD — Proportions (nLung vs tLung) ===\n")
print(results_luad_prop %>% arrange(p_adj))

cat("\n=== LUAD — ModuleScore (nLung vs tLung) ===\n")
print(luad_score_result)

cat("\n=== ImmunoT — Proportions (MPR vs NMPR) ===\n")
print(results_immuno_prop %>% arrange(p_adj))

cat("\n=== ImmunoT — ModuleScore (MPR vs NMPR) ===\n")
print(immuno_score_result)

# Save tables
write.csv(results_luad_prop,   file.path(output_dir, "22a_MixedEffects_LUAD_proportions.csv"),   row.names = FALSE)
write.csv(results_immuno_prop, file.path(output_dir, "22b_MixedEffects_ImmunoT_proportions.csv"), row.names = FALSE)

# ============================================================
# VISUALIZATION — Forest plot with CI
# ============================================================
results_all <- bind_rows(results_luad_prop, results_immuno_prop) %>%
  mutate(
    CI_low  = estimate - 1.96 * se,
    CI_high = estimate + 1.96 * se,
    sig     = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    dataset = factor(dataset, levels = c("LUAD", "ImmunoT"))
  )

p_forest <- ggplot(results_all,
                   aes(x = estimate, y = cluster,
                       color = dataset, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                 height = 0.3, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c("LUAD" = "#4575B4", "ImmunoT" = "#D73027")) +
  scale_shape_manual(values = c("***" = 16, "**" = 17, "*" = 15, "ns" = 1),
                     name = "Significance") +
  facet_wrap(~ dataset, scales = "free_x") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y      = element_text(size = 11, color = "#000000"),
    axis.text.x      = element_text(size = 11, color = "#000000"),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "#F5F5F5"),
    legend.position  = "right",
    plot.title       = element_text(size = 14, face = "bold")
  ) +
  labs(
    x     = "Estimate (effect of condition on proportion)",
    y     = "",
    title = "Mixed-effects model — CD8 cluster proportions",
    subtitle = "Patient as random effect | BH correction"
  )

ggsave(file.path(fig_dir, "22a_MixedEffects_ForestPlot.png"),
       p_forest, width = 14, height = 7, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "22a_MixedEffects_ForestPlot.pdf"),
       p_forest, width = 14, height = 7, bg = "white")

cat("\nDone — Script 22 complete\n")
