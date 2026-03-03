# CD8 T Cell Functional States in Lung Cancer — A Two-Dataset scRNA-seq Portfolio

**Author:** Myriam Yasmina Soumahoro  
**Master in Biology — University of Geneva**  
**Period:** January – February 2026  
**Environment:** R 4.4.1 (RStudio) — Seurat, Slingshot, CellChat  

---

## Biological Question

Do CD8 T cell functional states reshaping in the lung tumor microenvironment predict immunotherapy response in NSCLC?

This portfolio addresses this question through two independent scRNA-seq analyses, connected by a cross-dataset synthesis based on shared CD8 exhaustion signatures.

---

## Project Structure

```
CD8_NSCLC_scRNAseq/
├── README.md                    ← README global
├── LUAD_TME_CD8/
│   ├── README.md                ← README LUAD
│   ├── Data/
│   ├── Objects/
│   ├── Scripts/
│   └── Results/
│       ├── figures/
│       └── markers/
└── NSCLC_Immunotherapy_CD8/
    ├── README.md                ← README Dataset 2
    ├── Data/
    ├── Objects/
    ├── Scripts/
    └── Results/
        ├── figures/
        └── markers/
```

---

## Dataset 1 — LUAD Tumor Microenvironment (GSE131907)

**Publication:** Kim N. et al. Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. *Nature Communications* 2020;11:2285. https://doi.org/10.1038/s41467-020-16164-1  
**Platform:** 10x Genomics Chromium

### Dataset Composition

| Parameter | Value |
|---|---|
| Samples | 58 |
| Compartments | 7 (tLung, nLung, PE, nLN, mLN, tL/B, mBrain) |
| Total cells (post-QC) | ~41,000 (downsampled from >200,000) |
| Focus | nLung vs. tLung |

### CD8 T Cell States

Four states identified through unsupervised clustering, annotated based on top 50 differentially expressed markers:

| Cluster | Key Markers | nLung (n) | tLung (n) |
|---|---|---|---|
| CD8_Naive_CM | TCF7, IL7R, FOXP1, LTB | 245 | 243 |
| CD8_Effector_GZMK | GZMK, NKG7, CST7, TIGIT | 132 | 252 |
| CD8_TRM_Cytotoxic | GZMB, ITGAE, CXCR6, ZNF683, TIGIT | 17 | 223 |
| CD8_Proliferating | MKI67, TOP2A, CDK1 | 13 | 34 |

### Analysis Pipeline

```
Raw data (GSE131907)
        ↓
Quality Control (nFeature, nCount, percent.mt filtering)
        ↓
Normalization + Scaling (Seurat)
        ↓
PCA → Clustering (res = 0.5) → UMAP
        ↓
TME Annotation — 14 major cell types
        ↓
T cell Subset → 9 functional states annotated
        ↓
CD8 Subset → 4 states retained
        ↓
Trajectory Analysis (Slingshot)
        ↓
FFunctional Scoring (AddModuleScore)
— Exhaustion Terminal Score (PDCD1, TIGIT, HAVCR2, LAG3, CTLA4, ENTPD1, LAYN, PRDM1)
        ↓
Cell-Cell Communication (CellChat — nLung vs. tLung)
```

### Key Observations
- CD8_TRM_Cytotoxic cells appear markedly enriched in tLung vs. nLung (17 → 223 cells), though the low nLung count limits formal comparison
- GZMB and TIGIT expression tend to be higher in tLung
- Pseudotime inference suggests a possible differentiation continuum from Naive_CM toward TRM_Cytotoxic states
- CD8_TRM_Cytotoxic–to–Myeloid signaling is broader in tLung (1488 vs. 1103 interactions), notably via IFNG–IFNGR, CCL5–CCR1, HLA class II, PPIA–BSG, and MIF–CD74 axes

### Analytical Notes and Limitations
- Classical exhaustion markers (PDCD1, LAG3, HAVCR2, TOX) were not detectable at sufficient levels; narrative is reframed around TRM dysfunction (TIGIT+GZMB+)
- CD8_TRM_Cytotoxic had only 17 cells in nLung — comparisons are interpreted descriptively
- CellChat was run on nLung and tLung only, excluding other compartments

---

## Dataset 2 — Anti-PD-1 Immunotherapy Response in NSCLC (GSE207422)

**Publication:** Hu J, Zhang L, Xia H, Yan Y et al. Tumor microenvironment remodeling after neoadjuvant immunotherapy in non-small cell lung cancer revealed by single-cell RNA sequencing. *Genome Med* 2023;15:14. PMID: 36869384  
**Platform:** 10x Genomics Chromium

### Dataset Composition

| Condition | Samples | CD8+ T cells |
|---|---|---|
| PR (Major Pathologic Response) | 10 | 15,970 |
| SD (Non-Major Pathologic Response) | 6 | 4,349 |
| **Total** | **15** | **20,319** |

> 15 patients total (3 pre-treatment + 12 post-treatment). Analysis focuses on post-treatment samples stratified by pathologic response.

### CD8 T Cell States

Seven states identified through unsupervised clustering, annotated based on top 50 differentially expressed markers:

| Cluster | Function |
|---|---|
| CD8_Effector_GZMK | Cytotoxic effector |
| CD8_Exhausted_Terminal | Terminal exhaustion |
| CD8_Terminal_CX3CR1 | Terminal exhaustion — migratory |
| CD8_Proliferating | Active proliferation |
| CD8_TRM_like | Tissue-resident memory |
| CD8_Early_Activated_NR4A_high | Early activation |
| CD8_Activated_HLAII_high | MHC-II-activated |
| CD8_IFN_Stress_Response | Interferon stress response |

> For cross-dataset comparisons with LUAD, clusters were mapped to 4 equivalent states:
> CD8_Effector_GZMK, CD8_Exhausted_Terminal, CD8_TRM_like, CD8_Proliferating.

### Analysis Pipeline

```
Raw data (GSE207422)
        ↓
Quality Control + Normalization (Seurat)
        ↓
T cell Subclustering → CD8 Annotation (top 50 markers)
        ↓
Differential Abundance — PR vs. SD (barplot)
        ↓
Pseudotime Trajectory (Slingshot)
        ↓
Exhaustion Terminal Module Score
(PDCD1, TIGIT, HAVCR2, LAG3, CTLA4, ENTPD1, LAYN, PRDM1)
        ↓
Cell-Cell Communication (CellChat — PR vs. SD)
        ↓
Cross-Dataset Synthesis (TIGIT, PDCD1, Exhaustion Score)
```

### Key Observations
- CD8_Exhausted_Terminal cells appear more enriched in SD vs. PR samples
- Exhaustion Terminal Score tends to be higher in SD, suggesting a possible link between terminal exhaustion and reduced response to anti-PD-1
- CD8_Exhausted_Terminal cells show stronger predicted interactions with TAMs in SD (MHC-I, CCL5-CCR1 axes)
- These observations are exploratory and would require larger cohorts and functional validation
- CellChat analysis suggests slightly stronger overall interactions in SD vs. PR (269.2 vs. 250.8)

### Analytical Notes and Limitations
- Small and unbalanced patient numbers (PR n=10 samples, SD n=6 samples)
- No TCR clonotype or functional validation data available
- CD8 annotation relies on transcriptional markers only

---

## Cross-Dataset Synthesis

To bridge the two analyses, TIGIT and PDCD1 expression were visualized across all four conditions (nLung / tLung / PR / SD) in a unified 4-panel FeaturePlot. An Exhaustion Terminal Module Score was computed using the same gene signature in both datasets.

**Shared exhaustion signature genes:** PDCD1, TIGIT, HAVCR2, LAG3, CTLA4, ENTPD1, LAYN, PRDM1

The convergence of expression patterns between tLung (LUAD) and SD (immunotherapy) suggests that tumor-associated CD8 exhaustion observed in the LUAD microenvironment may carry functional relevance for immunotherapy response. This interpretation remains hypothesis-generating.

---

## Main Packages

| Package | Use |
|---|---|
| Seurat 5.x | Preprocessing, clustering, annotation |
| Slingshot (Bioconductor) | Pseudotime trajectory |
| CellChat (jinworks) | Cell-cell communication |
| ggplot2 / patchwork | Visualization |
