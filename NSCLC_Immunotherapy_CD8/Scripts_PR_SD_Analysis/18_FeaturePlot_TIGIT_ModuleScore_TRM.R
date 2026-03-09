#!/usr/bin/env Rscript
# ============================================================
# 18_FeaturePlot_TIGIT_ModuleScore.R
# Synthesis — TIGIT / PDCD1 expression + Exhaustion Terminal module score
# LUAD (nLung vs tLung) + Immunotherapy (PR vs SD)
# ============================================================
library(Seurat)
library(ggplot2)
library(patchwork)

fig_dir <- "Results/figures"

# -----------------------------
# Load objects
# -----------------------------
seu_luad   <- readRDS("Objects/15_CD8_annotated_states.rds")
seu_immuno <- readRDS("Objects/07_CD8_annotated.rds")

# Split immuno by response
seu_PR <- subset(seu_immuno, subset = RECIST == "PR")
seu_SD <- subset(seu_immuno, subset = RECIST == "SD")

# Split LUAD by tissue
seu_nLung <- subset(seu_luad, subset = Sample_Origin == "nLung")
seu_tLung <- subset(seu_luad, subset = Sample_Origin == "tLung")

# -----------------------------
# Exhaustion Terminal signature
# -----------------------------
exhaustion_genes <- list(Exhaustion_Terminal_signature = c("PDCD1", "TIGIT", "HAVCR2",
                                                           "LAG3", "CTLA4", "ENTPD1",
                                                           "LAYN", "PRDM1"))

# Compute module scores
seu_PR    <- AddModuleScore(seu_PR,    features = exhaustion_genes, name = "Exhaustion_score")
seu_SD    <- AddModuleScore(seu_SD,    features = exhaustion_genes, name = "Exhaustion_score")
seu_nLung <- AddModuleScore(seu_nLung, features = exhaustion_genes, name = "Exhaustion_score")
seu_tLung <- AddModuleScore(seu_tLung, features = exhaustion_genes, name = "Exhaustion_score")

# -----------------------------
# Plot 1 — FeaturePlot TIGIT 4 panels
# -----------------------------
p1 <- FeaturePlot(seu_nLung, features = "TIGIT", pt.size = 0.3,
                  min.cutoff = 0, max.cutoff = 3) +
  ggtitle("LUAD — Normal Lung") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p2 <- FeaturePlot(seu_tLung, features = "TIGIT", pt.size = 0.3,
                  min.cutoff = 0, max.cutoff = 3) +
  ggtitle("LUAD — Tumor Lung") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p3 <- FeaturePlot(seu_PR, features = "TIGIT", pt.size = 0.3,
                  min.cutoff = 0, max.cutoff = 3) +
  ggtitle("Immunotherapy — PR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p4 <- FeaturePlot(seu_SD, features = "TIGIT", pt.size = 0.3,
                  min.cutoff = 0, max.cutoff = 3) +
  ggtitle("Immunotherapy — SD") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_tigit <- (p1 | p2 | p3 | p4) +
  plot_annotation(title = "TIGIT expression across conditions",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(fig_dir, "18a_FeaturePlot_TIGIT_4panels.png"),
       p_tigit, width = 20, height = 5, dpi = 450, bg = "white")

# -----------------------------
# Plot 2 — FeaturePlot PDCD1 4 panels
# -----------------------------
p1b <- FeaturePlot(seu_nLung, features = "PDCD1", pt.size = 0.3,
                   min.cutoff = 0, max.cutoff = 3) +
  ggtitle("LUAD — Normal Lung") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p2b <- FeaturePlot(seu_tLung, features = "PDCD1", pt.size = 0.3,
                   min.cutoff = 0, max.cutoff = 3) +
  ggtitle("LUAD — Tumor Lung") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p3b <- FeaturePlot(seu_PR, features = "PDCD1", pt.size = 0.3,
                   min.cutoff = 0, max.cutoff = 3) +
  ggtitle("Immunotherapy — PR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p4b <- FeaturePlot(seu_SD, features = "PDCD1", pt.size = 0.3,
                   min.cutoff = 0, max.cutoff = 3) +
  ggtitle("Immunotherapy — SD") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_pdcd1 <- (p1b | p2b | p3b | p4b) +
  plot_annotation(title = "PDCD1 (PD-1) expression across conditions",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(fig_dir, "18b_FeaturePlot_PDCD1_4panels.png"),
       p_pdcd1, width = 20, height = 5, dpi = 450, bg = "white")

# -----------------------------
# Plot 3 — Exhaustion Terminal Score — PR vs SD
# -----------------------------
p_exh_PR <- FeaturePlot(seu_PR, features = "Exhaustion_score1", pt.size = 0.3,
                        cols = c("lightgrey", "#D73027"),
                        min.cutoff = -0.5, max.cutoff = 2) +
  ggtitle("Exhaustion Terminal Score — PR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_exh_SD <- FeaturePlot(seu_SD, features = "Exhaustion_score1", pt.size = 0.3,
                        cols = c("lightgrey", "#D73027"),
                        min.cutoff = -0.5, max.cutoff = 2) +
  ggtitle("Exhaustion Terminal Score — SD") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_exh_immuno <- p_exh_PR | p_exh_SD

ggsave(file.path(fig_dir, "18c_ModuleScore_Exhaustion_PRvsSD.png"),
       p_exh_immuno, width = 12, height = 5, dpi = 450, bg = "white")

# -----------------------------
# Plot 4 — Exhaustion Terminal Score — LUAD nLung vs tLung
# -----------------------------
p_exh_nLung <- FeaturePlot(seu_nLung, features = "Exhaustion_score1", pt.size = 0.3,
                           cols = c("lightgrey", "#D73027"),
                           min.cutoff = -0.5, max.cutoff = 1) +
  ggtitle("Exhaustion Terminal Score — Normal Lung") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_exh_tLung <- FeaturePlot(seu_tLung, features = "Exhaustion_score1", pt.size = 0.3,
                           cols = c("lightgrey", "#D73027"),
                           min.cutoff = -0.5, max.cutoff = 1) +
  ggtitle("Exhaustion Terminal Score — Tumor Lung") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_exh_luad <- p_exh_nLung | p_exh_tLung

ggsave(file.path(fig_dir, "18d_ModuleScore_Exhaustion_LUAD.png"),
       p_exh_luad, width = 12, height = 5, dpi = 450, bg = "white")

cat("\nScript 18 — all figures saved\n")