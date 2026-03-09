#!/usr/bin/env Rscript
# ============================================================
# 21_CellChat_CrossDataset_Comparison.R
# CellChat interaction comparison — LUAD vs ImmunoT
# Probabilités brutes — comparaison directionnelle
# Note: datasets indépendants, comparaison qualitative uniquement
# ============================================================

library(CellChat)
library(ggplot2)
library(dplyr)

fig_dir <- "Results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load CellChat objects
# -----------------------------
cc_nLung <- readRDS("../LUAD/objects/19_CellChat_nLung.rds")
cc_tLung <- readRDS("../LUAD/objects/19_CellChat_tLung.rds")
cc_PR    <- readRDS("../Immunotherapy/objects/12_CellChat_PR.rds")
cc_SD    <- readRDS("../Immunotherapy/objects/12_CellChat_SD.rds")

# -----------------------------
# Interactions of interest
# -----------------------------
interactions_of_interest <- c("CCL5 - CCR1", "PPIA - BSG")

extract_prob <- function(cc_obj, source_cluster, interactions, condition_label) {
  df <- subsetCommunication(cc_obj, sources.use = source_cluster)
  if (nrow(df) == 0) {
    message("No interactions found for: ", condition_label)
    return(NULL)
  }
  df %>%
    filter(interaction_name_2 %in% interactions) %>%
    group_by(interaction_name_2) %>%
    summarise(prob_mean = mean(prob, na.rm = TRUE)) %>%
    mutate(condition = condition_label)
}

# LUAD
df_nLung <- extract_prob(cc_nLung, "CD8_TRM_Cytotoxic", interactions_of_interest, "nLung")
df_tLung <- extract_prob(cc_tLung, "CD8_TRM_Cytotoxic", interactions_of_interest, "tLung")

# ImmunoT
df_PR <- extract_prob(cc_PR, "CD8_Cytotoxic_Exhausted", interactions_of_interest, "PR")
df_SD <- extract_prob(cc_SD, "CD8_Cytotoxic_Exhausted", interactions_of_interest, "SD")

# Combine
df_all <- bind_rows(df_nLung, df_tLung, df_PR, df_SD)
df_all$condition <- factor(df_all$condition, levels = c("nLung", "tLung", "PR", "SD"))
df_all$dataset   <- ifelse(df_all$condition %in% c("nLung", "tLung"),
                           "LUAD (GSE131907)",
                           "Immunotherapy (GSE207422)")

# -----------------------------
# Colors
# -----------------------------
condition_colors <- c(
  "nLung" = "#4575B4",
  "tLung" = "#D73027",
  "PR"    = "#4575B4",
  "SD"    = "#D73027"
)

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(df_all, aes(x = condition, y = prob_mean, fill = condition)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = condition_colors) +
  facet_grid(interaction_name_2 ~ dataset, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 12, color = "#000000"),
    axis.text.y      = element_text(size = 11, color = "#000000"),
    strip.text       = element_text(size = 12, face = "bold", color = "#000000"),
    strip.background = element_rect(fill = "#F5F5F5"),
    legend.position  = "none",
    plot.title       = element_text(size = 14, face = "bold", color = "#000000"),
    plot.subtitle    = element_text(size = 11, color = "#555555", face = "italic"),
    plot.caption     = element_text(size = 11, color = "#888888", face = "italic")
  ) +
  labs(
    x        = "",
    y        = "Mean interaction probability (CellChat)",
    title    = "Shared CD8 interaction patterns across datasets",
    subtitle = "CCL5–CCR1 preserved in responders; PPIA–BSG enriched in non-responders",
    caption  = "Note: datasets are independent; cross-dataset comparison is qualitative only."
  )

ggsave(file.path(fig_dir, "21_CellChat_CrossDataset_Comparison.png"),
       p, width = 10, height = 7, dpi = 450, bg = "white")
ggsave(file.path(fig_dir, "21_CellChat_CrossDataset_Comparison.pdf"),
       p, width = 10, height = 7, bg = "white")

cat("\nDone\n")