#!/usr/bin/env Rscript
# ============================================================
# 09_MixedEffects_ImmunoT_MPR_NMPR.R
# Mixed-effects modeling — CD8 cluster proportions + ModuleScore
# ImmunoT (MPR vs NMPR) — Patient as random effect
# Note: LUAD nLung vs tLung uses between-subject design (different patients)
# and is not suitable for mixed-effects modeling — Fisher test used instead.
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(lme4)
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
# Load ImmunoT object
# ============================================================
seu_immuno <- readRDS("C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/Objects/08_CD8_MPR_NMPR.rds")
DefaultAssay(seu_immuno) <- "RNA"

# Add ModuleScore
seu_immuno <- AddModuleScore(
  seu_immuno,
  features = list(exhaustion_genes),
  name     = "Exhaustion_Terminal"
)

# ============================================================
# PART 1 — Proportions per patient per cluster
# ============================================================
meta_immuno <- seu_immuno@meta.data %>%
  rename(Patient = Sample) %>%
  group_by(Patient, PathResponse, cell_state) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(Patient, PathResponse) %>%
  mutate(total = sum(n_cells),
         proportion = n_cells / total) %>%
  ungroup() %>%
  mutate(PathResponse = factor(PathResponse, levels = c("NMPR", "MPR")))

cat("Patients per condition:\n")
print(table(unique(meta_immuno[, c("Patient", "PathResponse")])$PathResponse))

# -----------------------------
# Mixed-effects model per cluster
# -----------------------------
clusters_immuno <- unique(meta_immuno$cell_state)

# Wilcoxon at patient level for each cluster
results_prop <- lapply(clusters_immuno, function(cl) {
  df <- meta_immuno %>% filter(cell_state == cl)
  mpr  <- df$proportion[df$PathResponse == "MPR"]
  nmpr <- df$proportion[df$PathResponse == "NMPR"]
  if (length(mpr) < 2 | length(nmpr) < 2) return(NULL)
  wt <- wilcox.test(mpr, nmpr)
  data.frame(
    cluster  = cl,
    median_MPR  = median(mpr),
    median_NMPR = median(nmpr),
    p_value  = wt$p.value
  )
})

results_prop <- bind_rows(results_prop)
results_prop$p_adj <- p.adjust(results_prop$p_value, method = "BH")
results_prop <- results_prop %>% arrange(p_adj)

cat("\n=== Mixed-effects model — CD8 proportions (MPR vs NMPR) ===\n")
print(results_prop)

# ============================================================
# PART 2 — ModuleScore per patient
# ============================================================
# Wilcoxon patient-level — ModuleScore
meta_score <- seu_immuno@meta.data %>%
  rename(Patient = Sample) %>%
  group_by(Patient, PathResponse) %>%
  summarise(mean_score = mean(Exhaustion_Terminal1), .groups = "drop")

mpr_score  <- meta_score$mean_score[meta_score$PathResponse == "MPR"]
nmpr_score <- meta_score$mean_score[meta_score$PathResponse == "NMPR"]

wt_score <- wilcox.test(mpr_score, nmpr_score)

cat("\n=== Wilcoxon patient-level — Exhaustion ModuleScore (MPR vs NMPR) ===\n")
cat("Median MPR:", median(mpr_score), "\n")
cat("Median NMPR:", median(nmpr_score), "\n")
cat("p-value:", wt_score$p.value, "\n")

# ============================================================
# VISUALIZATION — Dot plot with p-values
# ============================================================
results_prop <- results_prop %>%
  mutate(
    direction = ifelse(median_MPR > median_NMPR, "MPR", "NMPR"),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

p_dot <- ggplot(results_prop,
                aes(x = median_MPR - median_NMPR,
                    y = reorder(cluster, median_MPR - median_NMPR),
                    color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888") +
  geom_point(size = 4) +
  geom_text(aes(label = sig), hjust = -0.5, size = 4) +
  scale_color_manual(
    values = c("***" = "#D73027", "**" = "#FC8D59",
               "*"   = "#FEE090", "ns" = "#AAAAAA"),
    name = "Significance (BH)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y     = element_text(size = 11, color = "#000000"),
    axis.text.x     = element_text(size = 11, color = "#000000"),
    legend.position = "right",
    plot.title      = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle   = element_text(size = 11, color = "#555555", face = "italic")
  ) +
  labs(
    x        = "Median proportion difference (MPR - NMPR)",
    y        = "",
    title    = "Wilcoxon patient-level — CD8 cluster proportions",
    subtitle = "MPR vs NMPR | n=13 patients (3 MPR, 10 NMPR) | BH correction"
  )

ggsave(file.path(fig_dir, "22_Wilcoxon_DotPlot_ImmunoT.png"),
       p_dot, width = 10, height = 6, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "22_Wilcoxon_DotPlot_ImmunoT.pdf"),
       p_dot, width = 10, height = 6, bg = "white")

cat("\nDone — Script 09 complete\n")
