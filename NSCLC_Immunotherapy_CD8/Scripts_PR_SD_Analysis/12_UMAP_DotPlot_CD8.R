# ============================================================
# CD8_UMAP_DotPlot
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

Idents(seu_cd8) <- "cell_state"
dir.create("Results/figures", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Palette et markers
# ============================================================
cd8_colors <- c(
  "CD8_Effector_GZMK"             = "#92C5DE",  
  "CD8_Exhausted_Terminal"        = "#D73027",  
  "CD8_Terminal_CX3CR1"           = "#FDBB84",  
  "CD8_Proliferating"             = "#CBC9E2",  
  "CD8_TRM_like"                  = "#D6A5A5",  
  "CD8_Early_Activated_NR4A_high" = "#A1D99B",  
  "CD8_Activated_HLAII_high"      = "#D9D97A",  
  "CD8_IFN_Stress_Response"       = "#A8D8E8"  
)

cluster_order <- names(cd8_colors)

markers_cd8 <- c(
  "GZMK", "TCF7", "NKG7",
  "ENTPD1", "HAVCR2", "LAG3", "TIGIT", "CXCL13",
  "CX3CR1", "KLRF1", "FGFBP2",
  "MKI67", "TOP2A",
  "ZNF683", "CXCR6", "ITGAE", "GZMB",
  "NR4A1", "NR4A3",
  "VCAM1", "HLA-DRB5",
  "ISG15", "MX1"
)

# ============================================================
# UMAP global
# ============================================================
p_umap_global <- DimPlot(
  seu_cd8,
  reduction  = "umap",
  group.by   = "cell_state",
  cols       = cd8_colors,
  label      = TRUE,
  label.size = 2.5,
  repel      = TRUE,
  pt.size    = 0.3
) +
  ggtitle("CD8 T cell states") +
  theme_classic(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(5, "mm"),           # augmenté de 4 à 5
    legend.spacing.y = unit(3, "mm"),           # ajouté
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 8),
    strip.text       = element_text(size = 10, face = "bold")
  )

ggsave("Results/figures/CD8_UMAP_global.png",
       p_umap_global,
       width = 8, height = 5, dpi = 300, bg = "white")

# ============================================================
# UMAP split PR vs SD
# ============================================================
p_umap_split <- DimPlot(
  seu_cd8,
  reduction  = "umap",
  group.by   = "cell_state",
  split.by   = "RECIST",
  cols       = cd8_colors,
  label      = TRUE,
  label.size = 2,
  repel      = TRUE,
  pt.size    = 0.3
) +
  ggtitle("CD8 T cell states — PR vs SD") +
  theme_classic(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(5, "mm"),           # augmenté de 4 à 5
    legend.spacing.y = unit(3, "mm"),           # ajouté
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 8),
    strip.text       = element_text(size = 10, face = "bold")
  )

ggsave("Results/figures/CD8_UMAP_split_PR_SD.png",
       p_umap_split,
       width = 12, height = 5, dpi = 300, bg = "white")  

# ============================================================
# DotPlot PR et SD séparés
# ============================================================
seu_cd8_PR <- subset(seu_cd8, subset = RECIST == "PR")
seu_cd8_SD <- subset(seu_cd8, subset = RECIST == "SD")

Idents(seu_cd8_PR) <- "cell_state"
Idents(seu_cd8_SD) <- "cell_state"

p_dot_PR <- DotPlot(
  seu_cd8_PR,
  features = markers_cd8,
  cols     = c("lightgrey", "#08306b"),
  dot.scale = 5
) +
  coord_flip() +
  scale_y_discrete(limits = cluster_order) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),   # clusters
    axis.text.y     = element_text(size = 10),                          # markers
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(4, "mm"),
    legend.title    = element_text(size = 9),
    axis.title      = element_text(size = 10)
  ) +
  guides(
    size  = guide_legend(order = 1, title = "pct.exp",
                         override.aes = list(size = c(1, 3, 5, 7))),
    color = guide_colorbar(order = 2, title = "avg.exp.scaled",
                           barwidth = 1, barheight = 5)
  )
  ggtitle("PR")

p_dot_SD <- DotPlot(
  seu_cd8_SD,
  features = markers_cd8,
  cols     = c("lightgrey", "#08306b"),
  dot.scale = 5
) +
  coord_flip() +
  scale_y_discrete(limits = cluster_order) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y     = element_text(size = 10),
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(4, "mm"),
    legend.title    = element_text(size = 9),
    axis.title      = element_text(size = 10)
  ) +
  guides(
    size  = guide_legend(order = 1, title = "pct.exp",
                         override.aes = list(size = c(1, 3, 5, 7))),
    color = guide_colorbar(order = 2, title = "avg.exp.scaled",
                           barwidth = 1, barheight = 5)
  )
  ggtitle("SD")

ggsave("Results/figures/CD8_DotPlot_PR.png",
       p_dot_PR,
       width = 9, height = 8, dpi = 450, bg = "white")

ggsave("Results/figures/CD8_DotPlot_SD.png",
       p_dot_SD,
       width = 9, height = 8, dpi = 450, bg = "white")

cat("Done\n")

cat("Done\n")