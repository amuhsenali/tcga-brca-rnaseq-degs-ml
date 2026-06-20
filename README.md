# TCGA-BRCA RNA-seq Differential Expression and Machine Learning Analysis

**Ahmed Mohsin Ali¹**

¹Department of Computer Science, Jamia Millia Islamia, Jamia Nagar, New Delhi, India, 110025

---

## 1. Overview

End-to-end RNA-seq analysis of TCGA-BRCA data covering QC, normalization, differential expression, functional enrichment, and ML-based biomarker ranking. 2,384 significant DEGs identified; Random Forest, LASSO, and SVM models achieved >96% classification accuracy (AUC ≈ 1.0).

## 2. Dataset

| Detail | Value |
|---|---|
| Source | TCGA-BRCA (RNA-seq, STAR counts) |
| Tumor / Normal / Total samples | 100 / 113 / 213 |
| Genes (raw) | 60,660 |
| Genes tested (post-filter) | 25,571 |

## 3. Pipeline

Acquisition → QC (STAR two-pass, GRCh38) → Normalization (DESeq2 size factors + TPM) → DEG analysis (|log2FC| ≥ 2, p ≤ 0.05) → Enrichment (GO/KEGG/DO via clusterProfiler) → ML (RF, LASSO, SVM)

## 4. Key Results

- **2,384 significant DEGs** (1,135 up / 1,249 down)
- Top pathways: cell cycle, mitotic spindle organization, DNA replication, PI3K–Akt, p53 signaling
- Model performance:

| Model | Accuracy | Sensitivity | Specificity | AUC |
|---|---|---|---|---|
| Random Forest | 0.9683 | 0.9394 | 1.000 | 1.000 |
| LASSO | 1.0000 | 1.0000 | 1.000 | 1.000 |
| SVM | 0.9841 | 0.9697 | 1.000 | 1.000 |

- Top biomarker candidates: **VEGFD, MMP11, COL10A1, UBE2T, NEK2**

## 5. Repository Structure

```
tcga-brca-rnaseq-degs-ml/
├── scripts/
│   └── tcga_brca_rnaseq_analysis.R
├── figures/
├── tables/
├── sessionInfo.txt
├── LICENSE
└── README.md
```

## 6. Requirements

R 4.4.2 + key packages: `DESeq2`, `edgeR`, `limma`, `clusterProfiler`, `enrichplot`, `DOSE`, `org.Hs.eg.db`, `TCGAbiolinks`, `ggplot2`, `pheatmap`, `ComplexHeatmap`, `caret`, `randomForest`, `glmnet`, `pROC`, `tidyverse`, `biomaRt`, `pathview`.

Full session details: [`sessionInfo.txt`](./sessionInfo.txt)

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2","edgeR","limma","clusterProfiler","enrichplot",
                       "DOSE","org.Hs.eg.db","TCGAbiolinks","ComplexHeatmap",
                       "pathview","biomaRt"))
install.packages(c("ggplot2","pheatmap","caret","randomForest","glmnet",
                    "pROC","tidyverse"))
```

## 7. How to Run

```bash
git clone https://github.com/amuhsenali/tcga-brca-rnaseq-degs-ml.git
cd tcga-brca-rnaseq-degs-ml
Rscript scripts/tcga_brca_rnaseq_analysis.R
```

> Raw TCGA-BRCA data is not redistributed here; the script downloads it directly via `TCGAbiolinks` from the [GDC Data Portal](https://portal.gdc.cancer.gov/).

## 8. Citation

> Ali, A. M. (2026). *TCGA-BRCA RNA-seq Differential Expression and Machine Learning Analysis* [Computer software]. GitHub. https://github.com/amuhsenali/tcga-brca-rnaseq-degs-ml

## 9. License

MIT License — see [`LICENSE`](./LICENSE).

## 10. Contact

**Ahmed Mohsin Ali**
