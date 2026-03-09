# ============================================================
# 06_ModuleScore_Exhaustion_MPR_NMPR.R
# Exhaustion Terminal Module Score — MPR vs NMPR
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load object
seu_mpr <- readRDS("Objects/08_CD8_MPR_NMPR.rds")

# Split by response
seu_MPR  <- subset(seu_mpr, subset = PathResponse == "MPR")
seu_NMPR <- subset(seu_mpr, subset = PathResponse == "NMPR")

# Exhaustion Terminal signature
exhaustion_genes <- list(Exhaustion_Terminal_signature = c("PDCD1", "TIGIT", "HAVCR2",
                                                           "LAG3", "CTLA4", "ENTPD1",
                                                           "LAYN", "PRDM1"))

# Compute module scores
seu_MPR  <- AddModuleScore(seu_MPR,  features = exhaustion_genes, name = "Exhaustion_score")
seu_NMPR <- AddModuleScore(seu_NMPR, features = exhaustion_genes, name = "Exhaustion_score")

# FeaturePlot TIGIT
p1 <- FeaturePlot(seu_MPR, features = "TIGIT", pt.size = 0.3,
                  min.cutoff = 0, max.cutoff = 3) +
  ggtitle("Immunotherapy — MPR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p2 <- FeaturePlot(seu_NMPR, features = "TIGIT", pt.size = 0.3,
                  min.cutoff = 0, max.cutoff = 3) +
  ggtitle("Immunotherapy — NMPR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_tigit <- (p1 | p2) +
  plot_annotation(title = "TIGIT expression — MPR vs NMPR",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(fig_dir, "06a_FeaturePlot_TIGIT_MPR_NMPR.png"),
       p_tigit, width = 12, height = 5, dpi = 450, bg = "white")

# FeaturePlot PDCD1
p3 <- FeaturePlot(seu_MPR, features = "PDCD1", pt.size = 0.3,
                  min.cutoff = 0, max.cutoff = 3) +
  ggtitle("Immunotherapy — MPR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p4 <- FeaturePlot(seu_NMPR, features = "PDCD1", pt.size = 0.3,
                  min.cutoff = 0, max.cutoff = 3) +
  ggtitle("Immunotherapy — NMPR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_pdcd1 <- (p3 | p4) +
  plot_annotation(title = "PDCD1 (PD-1) expression — MPR vs NMPR",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(fig_dir, "06b_FeaturePlot_PDCD1_MPR_NMPR.png"),
       p_pdcd1, width = 12, height = 5, dpi = 450, bg = "white")

# Exhaustion Module Score
p5 <- FeaturePlot(seu_MPR, features = "Exhaustion_score1", pt.size = 0.3,
                  cols = c("lightgrey", "#D73027"),
                  min.cutoff = -0.5, max.cutoff = 2) +
  ggtitle("Exhaustion Terminal Score — MPR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p6 <- FeaturePlot(seu_NMPR, features = "Exhaustion_score1", pt.size = 0.3,
                  cols = c("lightgrey", "#D73027"),
                  min.cutoff = -0.5, max.cutoff = 2) +
  ggtitle("Exhaustion Terminal Score — NMPR") +
  theme(plot.title = element_text(size = 12, face = "bold", color = "black"))

p_exh <- (p5 | p6) +
  plot_annotation(title = "Exhaustion Terminal Module Score — MPR vs NMPR",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(fig_dir, "06c_ModuleScore_Exhaustion_MPR_NMPR.png"),
       p_exh, width = 12, height = 5, dpi = 450, bg = "white")

cat("Done — all figures saved\n")