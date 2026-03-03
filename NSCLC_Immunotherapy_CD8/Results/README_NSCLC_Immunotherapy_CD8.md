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
- CD8_Exhausted_Terminal cells appear more enriched in SD compared to PR samples, consistent with a higher exhaustion burden in non-responders
- The Exhaustion Terminal Score tends to be higher in SD, suggesting a possible link between terminal exhaustion and reduced response to anti-PD-1
- TIGIT and PDCD1 expression patterns in SD resemble those observed in LUAD tumor tissue, pointing to shared exhaustion features across contexts
- CellChat analysis suggests slightly stronger overall interactions in SD vs. PR (269.2 vs. 250.8), with CD8_Exhausted_Terminal cells showing broader predicted outgoing signaling toward TAMs and monocytes in SD, notably via MHC-I and CCL5-CCR1 axes
-These observations are exploratory and would require larger cohorts and functional validation to draw firm conclusions

---

## Project Structure
```
Dataset2_Immunotherapy/

├── Scripts/
│   └── R scripts (numbered pipeline)
└── Results/
    ├── figures/
    └── markers/
```
