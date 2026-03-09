# Dataset 2 — Anti-PD-1 Immunotherapy Response in NSCLC (GSE207422)
**Period:** February – March 2026  

## Overview
Single-cell RNA sequencing analysis of the tumor microenvironment in non-small cell lung 
cancer (NSCLC) patients treated with neoadjuvant PD-1 blockade combined with chemotherapy.
Patients were stratified by pathologic response assessed histologically on the surgical 
specimen: Major Pathologic Response (MPR, ≤10% residual viable tumor cells) vs 
Non-Major Pathologic Response (NMPR, >10% residual viable tumor cells), as defined 
in Hu et al. 2023.

> Initial analysis used RECIST radiological classification (PR/SD). This variable was 
> found unsuitable for a neoadjuvant context — pathologic response provides a more 
> biologically meaningful stratification. All analyses were repeated using MPR/NMPR.

**Publication:** Hu J, Zhang L, Xia H, Yan Y et al. Tumor microenvironment remodeling after 
neoadjuvant immunotherapy in non-small cell lung cancer revealed by single-cell RNA sequencing.
*Genome Med* 2023;15:14. PMID: 36869384  
**GEO Accession:** GSE207422

---

## Dataset Composition
| Condition | Patients | CD8+ T cells |
|-----------|----------|--------------|
| MPR       | 3        | 5,493        |
| NMPR      | 10       | 12,983       |
| **Total** | **13**   | **18,476**   |

> 15 patients total (3 pre-treatment excluded + 2 excluded: 1 non-evaluable, 1 pCR).  
> This analysis focuses on 13 post-treatment evaluable samples stratified by pathologic response.

---

## CD8+ T Cell Clusters
Eight functional states were identified through unsupervised clustering and annotated 
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

> For cross-dataset comparisons with LUAD, four equivalent states were retained:
> CD8_Effector_GZMK, CD8_Exhausted_Terminal, CD8_TRM_like, CD8_Proliferating.

---

## Analysis Pipeline

### 1. Quality Control & Preprocessing
- Doublet removal, mitochondrial content filtering
- Normalization, highly variable gene selection, PCA, UMAP

### 2. Pathologic Response Metadata Mapping
- MPR/NMPR labels mapped to CD8 object (`08_CD8_MPR_NMPR.rds`)
- MPR/NMPR labels mapped to TME object (`04_TME_MPR_NMPR.rds`)

### 3. T Cell Subclustering & CD8 Annotation
- Subclustering of T cells and identification of CD8+ populations
- Annotation based on top 50 markers per cluster
- Validation via DotPlot and FeaturePlot

### 4. UMAP Visualization
- UMAP colored by cluster and by condition (MPR/NMPR)
- UMAP split by pathologic response with cluster annotations

### 5. Differential Abundance — MPR vs NMPR
- Barplot of CD8 cluster proportions by pathologic response
- Fisher's exact test with BH correction on cluster proportions
- Wilcoxon patient-level test on cluster proportions (n=13 patients)
  — no significant differences detected, likely reflecting low statistical power (MPR n=3)

### 6. Module Scores & FeaturePlots
- Exhaustion Terminal Score (PDCD1, TIGIT, HAVCR2, LAG3, CTLA4, ENTPD1, LAYN, PRDM1)
- FeaturePlot PDCD1 and TIGIT — MPR vs NMPR
- Global ModuleScore: no significant difference MPR vs NMPR (Wilcoxon p=0.73),
  consistent with exhaustion being a sub-population phenomenon

### 7. Cell-Cell Communication — CellChat
- Separate CellChat objects constructed for MPR and NMPR (`04_TME_MPR_NMPR.rds`)
- Comparison of interaction strength and signaling patterns
- Bubble plot: CD8_Cytotoxic_Exhausted outgoing interactions (prob > 0.05)
- Note: interaction metrics not normalized by cell count given structural 
  group imbalance (MPR: ~5,493 cells, NMPR: ~12,983 cells)

### 8. Cross-Dataset Synthesis
- FeaturePlot PDCD1 across 4 conditions (nLung / tLung / MPR / NMPR)
- Exhaustion Terminal Module Score compared across LUAD and ImmunoT contexts
  (Wilcoxon cross-dataset: LUAD p<0.001, ImmunoT p=0.73)
- Schematic bifurcating model: CD8 remodeling under anti-PD-1 treatment

---

## Key Findings
- CD8_Exhausted_Terminal cells are significantly enriched in MPR (Fisher's exact test, 
  OR=3.36, p_adj<0.001), suggesting that exhausted-phenotype cells capable of PD-1 
  reactivation are more prevalent in pathologic responders
- CD8_Effector_GZMK also enriched in MPR (OR=1.55, p_adj<0.001), indicating a 
  co-occurrence of effector and exhausted states in responders
- CD8_IFN_Stress_Response, CD8_Early_Activated_NR4A_high, CD8_Activated_HLAII_high, 
  and CD8_Proliferating are significantly enriched in NMPR (p_adj<0.001), pointing to 
  dysfunctional and activation-impaired states in non-responders
- The Exhaustion Terminal Score does not show a significant global difference between 
  MPR and NMPR (p=0.73) — exhaustion is a sub-population phenomenon captured at the 
  cluster level but diluted across all CD8 states
- CellChat: stronger predicted interactions between CD8_Exhausted_Terminal and immunosuppressive TAMs in NMPR via PTPRC-MRC1; PPIA-BSG interactions shift from TAM_like in MPR to Tumor_epithelial in NMPR; CCL5-CCR1 preserved in both conditions 
  cells and immunosuppressive TAMs in NMPR, notably via PPIA-BSG and PTPRC-MRC1 axes
- CCL5-CCR1 interactions are preserved in both conditions

---

## Limitations
- Small patient cohort (n=13 evaluable), with marked group imbalance (MPR n=3, NMPR n=10)
- Patient-level statistical tests (Wilcoxon) underpowered due to low MPR sample size
- Cell-cell interaction metrics not normalized by cell count given structural imbalance
- CD8_Exhausted_Terminal annotation based on transcriptional markers only — 
  not validated at protein level or by functional assays
- Cross-dataset comparisons are qualitative and hypothesis-generating

---

## Project Structure
```
Immunotherapy/
├── Data/
│   └── raw count matrices (GSE207422)
├── Objects/
│   └── Seurat objects (.rds)
├── Scripts/
│   ├── 01_Add_PathResponse_metadata.R
│   ├── 01b_Add_PathResponse_TME.R
│   ├── 02_UMAP_MPR_NMPR.R
│   ├── 02b_UMAP_TME_MPR_NMPR.R
│   ├── 03_Barplot_TME_MPR_NMPR.R
│   ├── 04_Barplot_CD8_MPR_NMPR.R
│   ├── 05_ProportionTest_CD8_MPR_NMPR.R
│   ├── 06_ModuleScore_Exhaustion_MPR_NMPR.R
│   ├── 07_CellChat_MPR_NMPR_Create.R
│   └── 07b_CellChat_MPR_NMPR_Visualization.R
└── Results/
    ├── figures/
    └── tables/
```

## Main Packages
| Package | Version | Use |
|---------|---------|-----|
| Seurat | 5.x | Preprocessing, clustering, annotation |
| CellChat | jinworks/CellChat | Cell-cell communication |
| lme4 / lmerTest | CRAN | Mixed-effects modeling (attempted, underpowered) |
| ggplot2 / patchwork | CRAN | Visualization |

## Reference
Hu J. et al. (2023). Tumor microenvironment remodeling after neoadjuvant immunotherapy 
in non-small cell lung cancer revealed by single-cell RNA sequencing. 
*Genome Medicine*, 15, 14. https://doi.org/10.1186/s13073-023-01164-9
