# tcga-brca-rnaseq-degs-ml
End-to-end RNA-seq pipeline on TCGA-BRCA (213 samples): DESeq2 differential expression, clusterProfiler functional enrichment, and ML-based biomarker ranking (Random Forest, LASSO, SVM) achieving >96% classification accuracy

Requirements
R 4.4.2 + key packages: DESeq2, edgeR, limma, clusterProfiler, enrichplot, DOSE, org.Hs.eg.db, TCGAbiolinks, ggplot2, pheatmap, ComplexHeatmap, caret, randomForest, glmnet, pROC, tidyverse, biomaRt, pathview.
