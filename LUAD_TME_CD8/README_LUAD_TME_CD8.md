# scRNA-seq Analysis of the Tumor Microenvironment in Lung Adenocarcinoma
## Focus: CD8 T Cell Remodeling in Tumor vs. Normal Lung

**Author:** Myriam Yasmina Soumahoro  
**Dataset:** GSE131907 — Kim et al., 2020 (*Nature Communications*)  
**Period:** January – February 2026  
**Platform:** 10x Genomics Chromium  
**Environment:** R 4.4.1 (RStudio) — Seurat, Slingshot, CellChat  

---

## Biological Question

What functional states do CD8 T cells adopt in the lung tumor microenvironment (tLung) compared to normal lung tissue (nLung), and how does this remodeling shape intercellular communication?

---

## Dataset Overview

| Parameter | Value |
|---|---|
| Source | GEO: GSE131907 |
| Samples | 58 samples |
| Compartments | 7 (tLung, nLung, PE, nLN, mLN, tL/B, mBrain) |
| Total cells (post-QC) | ~41,000 (downsampled from >200,000) |
| Technology | 10x Genomics Chromium (scRNA-seq) |
| Focus | nLung vs. tLung |

---

## Analysis Workflow

```
Raw data (GSE131907)
        ↓
Quality Control
(nFeature, nCount, percent.mt filtering)
        ↓
Normalization + Scaling (Seurat)
        ↓
PCA → Clustering (res = 0.5) → UMAP
        ↓
TME Annotation — 14 major cell types
(7 compartments profiled)
        ↓
T cell Subset → 9 functional states annotated
(CD4, CD8, Treg, NK-like, Proliferating)
        ↓
CD8 Subset → 4 states retained
(Naive_CM, Effector_GZMK, TRM_Cytotoxic, Proliferating)
        ↓
Trajectory Analysis (Slingshot)
        ↓
Functional Scoring (AddModuleScore)
— Exhaustion Terminal Score (PDCD1, TIGIT, HAVCR2, LAG3, CTLA4, ENTPD1, LAYN, PRDM1)
        ↓
Cell-Cell Communication (CellChat — nLung vs. tLung)
```

---

## CD8 T Cell States

| Cluster | Annotation | Key Markers | nLung (n) | tLung (n) |
|---|---|---|---|---|
| 0 | CD8_Naive_CM | TCF7, IL7R, FOXP1, LTB | 245 | 243 |
| 1 | CD8_Effector_GZMK | GZMK, NKG7, CST7, TIGIT | 132 | 252 |
| 2 | CD8_TRM_Cytotoxic | GZMB, ITGAE, CXCR6, ZNF683, TIGIT | 17 | 223 |
| 4 | CD8_Proliferating | MKI67, TOP2A, CDK1 | 13 | 34 |

---

## Key Findings

**TME remodeling appears associated with a functional shift in CD8 T cells from nLung to tLung.**

These observations suggest a possible tumor-associated enrichment of TRM_Cytotoxic cells (17 → 223), though the low cell count in nLung limits formal comparison. GZMB and TIGIT expression appear elevated in tLung, consistent with a dysfunctional tissue-resident profile. Pseudotime analysis is compatible with a differentiation continuum toward TRM_Cytotoxic states, and the Exhaustion Terminal Score tends to be higher in tLung. CellChat analysis indicates broader CD8_TRM_Cytotoxic–Myeloid signaling in tLung (1488 vs. 1103 interactions), notably via IFNG–IFNGR, CCL5–CCR1, HLA class II, PPIA–BSG, and MIF–CD74 axes. These findings are descriptive and would require functional validation to draw firm conclusions.

---

## Analytical Notes and Limitations

- Classical exhaustion markers (PDCD1, LAG3, HAVCR2, TOX) were not detectable at sufficient levels in this dataset. The narrative is reframed around TRM dysfunction (TIGIT+GZMB+) rather than canonical exhaustion.
- The CD8_TRM_Cytotoxic cluster had only 17 cells in nLung, limiting the statistical power of nLung vs. tLung comparisons for this cluster. Observations are interpreted descriptively.
- CD8_TRM and CD8_TRM_CXCL13, identified as distinct states in the full T cell analysis, were merged into CD8_TRM_Cytotoxic in the CD8-focused analysis based on shared canonical markers (ITGAE, CXCR6, GZMB).
- CellChat was run on nLung and tLung compartments only (n ≈ 11,000 cells), excluding the 5 other compartments profiled.
- The two datasets (LUAD TME and Immunotherapy) are independent — the cross-dataset bridge is based on signature projection and should be interpreted as hypothesis-generating.

---

## Repository Structure

```

LUAD_TME_CD8/
├── README.md
├── Scripts/
├── └── Results/
│       ├── figures/
│       └── markers/

```

---

## Main Packages

| Package | Version | Use |
|---|---|---|
| Seurat | 5.x | Preprocessing, clustering, annotation |
| Slingshot | Bioconductor | Pseudotime trajectory |
| CellChat | jinworks/CellChat | Cell-cell communication |
| ggplot2 / patchwork | CRAN | Visualization |
| ggpubr | CRAN | Statistical annotation |

---

## Reference

Kim N. et al. (2020). Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. *Nature Communications*, 11, 2285. https://doi.org/10.1038/s41467-020-16164-1
