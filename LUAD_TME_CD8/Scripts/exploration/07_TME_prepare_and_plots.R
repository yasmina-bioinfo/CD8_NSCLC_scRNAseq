#!/usr/bin/env Rscript

# ============================================================
# 07_TME_prepare_and_plots.R
# LUAD (GSE131907) — Inject metadata + ordered factors + palette
# Outputs: UMAP major_celltype + TME composition barplots (prop + counts)
# Saves: PNG + PDF (white background)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# -----------------------------
# 0) Paths (EDIT THESE 2 LINES)
# -----------------------------
seurat_rds <- "objects/05_seurat_with_umap_clusters.rds"   # <-- your object
anno_tsv   <- "data/GSE131907_Lung_Cancer_cell_annotation.txt"

out_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Load Seurat object + annotation
# -----------------------------
obj  <- readRDS(seurat_rds)
anno <- read.delim(anno_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

stopifnot(is(obj, "Seurat"))
stopifnot(all(c("Barcode","Sample_Origin") %in% colnames(anno)))

# -----------------------------
# 2) Build a robust mapping: Seurat cellnames -> core barcode -> anno$Barcode
#    Seurat cellnames look like: "GGCG..._EFFUSION_13"
#    Core barcode is before the first "_"
# -----------------------------
cell_full <- colnames(obj)
cell_core <- sub("_.*$", "", cell_full)

# create mapping vectors
origin_map <- setNames(anno$Sample_Origin, anno$Barcode)

# inject Sample_Origin (as character first)
sample_origin <- origin_map[cell_core]

# sanity checks
na_n <- sum(is.na(sample_origin))
message("Cells in object: ", length(cell_full))
message("NAs after mapping Sample_Origin: ", na_n)

if (na_n == length(sample_origin)) {
  stop("Sample_Origin mapping failed (all NA). Check that barcodes match (core vs full) and anno_tsv path.")
}

# write metadata safely (NO overlap error, because rownames match EXACTLY)
meta_add <- data.frame(
  Sample_Origin = sample_origin,
  row.names = cell_full,
  stringsAsFactors = FALSE
)
obj <- AddMetaData(obj, metadata = meta_add)

# quick view
message("\nSample_Origin distribution (after injection):")
print(table(obj@meta.data$Sample_Origin, useNA = "ifany"))

# -----------------------------
# 3) Inject major cell type annotation
# -----------------------------

# Create named vectors from the annotation file
# (Barcode → Sample_Origin / Cell_type / Cell_type.refined)
origin_map   <- setNames(anno$Sample_Origin, anno$Barcode)
celltype_map <- setNames(anno$Cell_type, anno$Barcode)
refined_map  <- setNames(anno$Cell_type.refined, anno$Barcode)

# Extract core barcodes from Seurat object cell names
# (remove suffix after the first "_")
seu_cells <- colnames(obj)
seu_core  <- sub("_.*$", "", seu_cells)

# Inject metadata into the Seurat object
obj$Sample_Origin     <- unname(origin_map[seu_core])
obj$Cell_type         <- unname(celltype_map[seu_core])
obj$Cell_type_refined <- unname(refined_map[seu_core])

# Define a robust major_celltype variable:
# use Cell_type.refined when available,
# otherwise fall back to Cell_type
obj$major_celltype <- ifelse(
  is.na(obj$Cell_type_refined) | obj$Cell_type_refined == "",
  obj$Cell_type,
  obj$Cell_type_refined
)

# Quick sanity checks
cat("NAs Sample_Origin:", sum(is.na(obj$Sample_Origin)), "\n")
cat("NAs Cell_type:", sum(is.na(obj$Cell_type)), "\n")
cat("NAs major_celltype:", sum(is.na(obj$major_celltype)), "\n")

print(table(obj$major_celltype))

# -----------------------------
# 4) Define BIOLOGICAL order for compartments (Sample_Origin)
#    Keep only levels that exist in the object.
# -----------------------------
comp_order_full <- c("nLung", "tLung", "PE", "nLN", "mLN", "tL/B", "mBrain")
comp_present <- intersect(comp_order_full, unique(obj$Sample_Origin))
# (If there are any unexpected compartments, append them at the end)
comp_extra <- setdiff(unique(obj$Sample_Origin), comp_present)
comp_order <- c(comp_present, sort(comp_extra))

obj$Sample_Origin <- factor(obj$Sample_Origin, levels = comp_order)

# -----------------------------
# 5) Define STRATEGIC order for major cell types (stack + legend)
# -----------------------------
celltype_order_full <- c(
  "CD8 T cells",
  "CD4 T cells",
  "B lineage",
  "Macrophages",
  "Monocytes",
  "Dendritic cells",
  "Mast cells",
  "Fibroblasts",
  "Endothelial cells",
  "Lung epithelial",
  "Tumor basal/proliferative",
  "Tumor hypoxic/remodeling",
  "Tumor secretory/mucinous",
  "Other"
)

celltype_present <- intersect(celltype_order_full, unique(obj$major_celltype))
celltype_extra   <- setdiff(unique(obj$major_celltype), celltype_present)
celltype_order   <- c(celltype_present, sort(celltype_extra))

obj$major_celltype <- factor(obj$major_celltype, levels = celltype_order)

# -----------------------------
# 6) Palette (EDIT HERE if you want your own exact colors)
#    Tumor states intentionally more vivid.
# -----------------------------
# --- Define a COMPLETE named palette (names must match major_celltype EXACTLY) ---
tme_colors <- c(
  "CD8 T cells"                = "#1f77b4",
  "CD4 T cells"                = "#6baed6",
  "B lineage"                  = "#2ca02c",
  "Macrophages"                = "#ff7f0e",
  "Monocytes"                  = "#fdb462",
  "Dendritic cells"            = "#fb9a99",
  "Mast cells"                 = "#fdd0a2",
  "Fibroblasts"                = "#6a3d9a",
  "Endothelial cells"          = "#b15928",
  "Lung epithelial"            = "#9ecae1",
  "Tumor basal/proliferative"  = "#7f0000",  # vivid
  "Tumor hypoxic/remodeling"   = "#b30000",  # vivid
  "Tumor secretory/mucinous"   = "#e31a1c",  # vivid
  "Other"                      = "#bdbdbd"
)

# --- Safety checks (THIS is what prevents grey disasters) ---
ct_levels <- levels(obj$major_celltype)
missing_in_palette <- setdiff(ct_levels, names(tme_colors))
extra_in_palette   <- setdiff(names(tme_colors), ct_levels)

if (length(missing_in_palette) > 0) {
  stop("Missing colors for: ", paste(missing_in_palette, collapse = ", "))
}
if (length(extra_in_palette) > 0) {
  message("Note: unused colors in palette: ", paste(extra_in_palette, collapse = ", "))
}

# -----------------------------
# 7) UMAP — Major cell types (single UMAP you keep)
# -----------------------------
p_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "major_celltype",
  label     = FALSE,
  raster    = TRUE
) +
  scale_color_manual(values = tme_colors, drop = FALSE) +
  ggtitle("Tumor Microenvironment — Major Cell Types") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    legend.title = element_blank()
  )

ggsave(file.path(out_dir, "UMAP_major_celltype.png"), p_umap, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "UMAP_major_celltype.pdf"), p_umap, width = 10, height = 7, bg = "white")

# -----------------------------
# 8) Barplot — TME composition across compartments (PROPORTIONS)
# -----------------------------
df_comp <- obj@meta.data %>%
  filter(!is.na(Sample_Origin), !is.na(major_celltype)) %>%
  count(Sample_Origin, major_celltype, name = "n") %>%
  group_by(Sample_Origin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar_prop <- ggplot(df_comp, aes(x = Sample_Origin, y = prop, fill = major_celltype)) +
  geom_col(width = 0.85, color = NA) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = tme_colors, drop = FALSE, na.value = "grey80") +
  labs(
    title = "TME composition across compartments",
    x = "Compartment",
    y = "Proportion of cells"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(file.path(out_dir, "Barplot_TME_composition_prop.png"), p_bar_prop, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Barplot_TME_composition_prop.pdf"), p_bar_prop, width = 12, height = 8, bg = "white")

# -----------------------------
# 9) Barplot — Counts (optional but good QC)
# -----------------------------
p_bar_counts <- ggplot(df_comp, aes(x = Sample_Origin, y = n, fill = major_celltype)) +
  geom_col(width = 0.85, color = NA) +
  scale_fill_manual(values = tme_colors, drop = FALSE, na.value = "grey80") +
  labs(
    title = "TME composition across compartments (counts)",
    x = "Compartment",
    y = "Number of cells"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(file.path(out_dir, "Barplot_TME_composition_counts.png"), p_bar_counts, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Barplot_TME_composition_counts.pdf"), p_bar_counts, width = 12, height = 8, bg = "white")

# -----------------------------
# 10) Save updated object (now reproducible)
# -----------------------------
saveRDS(obj, file = "objects/06_seurat_TME_with_SampleOrigin_and_orderedCelltypes.rds")

message("\nDONE ✅ Figures saved in: ", out_dir)
message("Saved updated Seurat object: objects/06_seurat_TME_with_SampleOrigin_and_orderedCelltypes.rds")