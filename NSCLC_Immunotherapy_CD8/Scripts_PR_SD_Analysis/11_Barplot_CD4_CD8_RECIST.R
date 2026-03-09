# ============================================================
# Tcells_barplot.R
# Barplot CD4 and CD8 T cell states — PR vs SD
# Input : seu_cd4, seu_cd8 (objects already in memory)
# Output: Results/figures/*.png + *.pdf
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ============================================================
# Palettes
# ============================================================
cd4_colors <- c(
  "CD4_Naive_CentralMemory"       = "#6baed6",
  "CD4_Tfh_like"                  = "#2171b5",
  "Treg_Activated"                = "#f768a1",
  "CD4_Early_Activated_NR4A_high" = "#9ecae1",
  "Treg_Effector"                 = "#c51b8a"
)

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

# ============================================================
# Output directory
# ============================================================
outdir <- "Results/figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load objects
# ============================================================
seu_cd4 <- readRDS("Objects/07_CD4_annotated.rds")
seu_cd8 <- readRDS("Objects/07_CD8_annotated.rds")

# ============================================================
# Barplot theme — shared
# ============================================================
barplot_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.text       = element_text(size = 11),
    axis.title.y    = element_text(size = 11),
    legend.text     = element_text(size = 12),
    legend.key.size = unit(6, "mm"),
    legend.spacing.y = unit(0.5, "cm")
  )

# ============================================================
# Barplot function
# ============================================================
make_barplot <- function(seu_obj, color_pal, title) {
  
  meta <- seu_obj@meta.data
  
  df <- meta %>%
    filter(!is.na(cell_state), !is.na(RECIST)) %>%
    group_by(RECIST, cell_state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(RECIST) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()
  
  df$cell_state <- factor(df$cell_state, levels = names(color_pal))
  df$RECIST     <- factor(df$RECIST, levels = c("PR", "SD"))
  
  ggplot(df, aes(x = RECIST, y = pct, fill = cell_state)) +
    geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = color_pal, drop = FALSE) +
    scale_y_continuous(labels = percent_format(scale = 1),
                       expand = c(0, 0)) +
    labs(title = title, x = NULL, y = "Percentage of cells", fill = NULL) +
    barplot_theme
}

# ============================================================
# Generate plots
# ============================================================
p_cd4      <- make_barplot(seu_cd4, cd4_colors, "CD4 T cell states")
p_cd8      <- make_barplot(seu_cd8, cd8_colors, "CD8 T cell states")
p_combined <- p_cd4 | p_cd8

# ============================================================
# Save
# ============================================================
# CD8 only 
ggsave(file.path(outdir, "Tcells_barplot_CD8_PR_SD.png"),
       p_cd8, width = 8, height = 6, dpi = 450, bg = "white")
ggsave(file.path(outdir, "Tcells_barplot_CD8_PR_SD.pdf"),
       p_cd8, width = 8, height = 6, bg = "white")

# Combined CD4 + CD8 
ggsave(file.path(outdir, "Tcells_barplot_CD4_CD8_PR_SD.png"),
       p_combined, width = 14, height = 6, dpi = 450, bg = "white")
ggsave(file.path(outdir, "Tcells_barplot_CD4_CD8_PR_SD.pdf"),
       p_combined, width = 14, height = 6, bg = "white")

message("Done — barplots saved in: ", outdir)