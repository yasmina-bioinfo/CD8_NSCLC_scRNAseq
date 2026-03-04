# Dataset 2 — Anti-PD-1 Immunotherapy Response in NSCLC (GSE207422)

## Overview
Single-cell RNA sequencing analysis of the tumor microenvironment in non-small cell lung 
cancer (NSCLC) patients treated with neoadjuvant PD-1 blockade combined with chemotherapy.
Patients were stratified by pathological response: Major Pathologic Response (MPR/PR) 
vs Non-Major Pathologic Response (NMPR/SD).

**Publication:** Hu J, Zhang L, Xia H, Yan Y et al. Tumor microenvironment remodeling after 
neoadjuvant immunotherapy in non-small cell lung cancer revealed by single-cell RNA sequencing.
*Genome Med* 2023;15:14. PMID: 36869384  
**GEO Accession:** GSE207422

---

## Dataset Composition
| Condition | Samples | CD8+ T cells |
|-----------|---------|--------------|
| PR (MPR)  | 10      | 15,970       |
| SD (NMPR) | 6       | 4,349        |
| **Total** | **15**  | **20,319**   |

> 15 patients total (3 pre-treatment + 12 post-treatment).  
> This analysis focuses on post-treatment samples, stratified by pathologic response.

---

## CD8+ T Cell Clusters
Seven functional states were identified through unsupervised clustering and annotated 
based on the top 50 differentially expressed markers per cluster:

| Cluster | Function |
|---------|----------|
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

---

## Analysis Pipeline

### 1. Quality Control & Preprocessing
- Doublet removal, mitochondrial content filtering
- Normalization, highly variable gene selection, PCA, UMAP

### 2. T Cell Subclustering & CD8 Annotation
- Subclustering of T cells and identification of CD8+ populations
- Annotation based on top 50 markers per cluster
- Validation via DotPlot and FeaturePlot

### 3. Differential Abundance — PR vs SD
- Barplot of CD8 cluster proportions by pathologic response
- Fisher's exact test with BH correction on cluster proportions (PR vs SD)

### 4. Pseudotime Trajectory — Slingshot
- Lineage inference across CD8 functional states

### 5. Module Scores
- Exhaustion Terminal Score (PDCD1, TIGIT, HAVCR2, LAG3, CTLA4, ENTPD1, LAYN, PRDM1)

### 6. Cell-Cell Communication — CellChat
- Separate CellChat objects constructed for PR and SD
- Comparison of interaction strength and signaling patterns
- Bubble plot: CD8_Exhausted_Terminal outgoing interactions (prob > 0.05)

### 7. Cross-Dataset Synthesis
- FeaturePlot TIGIT and PDCD1 across 4 conditions (nLung / tLung / PR / SD)
- Exhaustion Terminal Module Score compared across LUAD and immunotherapy contexts

---

## Key findings
- CD8_Exhausted_Terminal cells appear more enriched in PR compared to SD samples 
  (Fisher's exact test, OR=32.1, p_adj<0.001), suggesting that exhausted-phenotype 
  cells capable of PD-1 reactivation are more prevalent in responders
- CD8_IFN_Stress_Response and CD8_Terminal_CX3CR1 are significantly enriched in SD 
  (p_adj<0.001), pointing to dysfunctional, stress-associated states in non-responders
- The Exhaustion Terminal Score does not show a consistent directional pattern between 
  PR and SD when computed across all CD8 cells, likely reflecting the higher abundance 
  of CD8_Exhausted_Terminal cells in PR
- CellChat analysis suggests stronger predicted interactions between CD8_Exhausted_Terminal 
  cells and immunosuppressive TAMs in SD, notably via PPIA-BSG and PTPRC-MRC1 axes
- CCL5-CCR1 interactions are preserved in PR, consistent with maintained effector signaling 
  in responders
- These observations are exploratory and would require larger cohorts and functional 
  validation to draw firm conclusions

---

## Project Structure
```
Dataset2_Immunotherapy/
├── Data/
│   └── raw count matrices (GSE207422)
├── Objects/
│   └── Seurat objects (.rds)
├── Scripts/
│   └── R scripts (numbered pipeline)
└── Results/
    ├── figures/
    └── markers/
```
