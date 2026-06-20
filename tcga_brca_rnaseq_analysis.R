# ============================================================================
# TCGA BRCA RNA-seq Differential Expression Analysis (Clean Pipeline)
# Optimized for MacBook Air M1 (8GB RAM)
# Course: NGS Data Analysis | Instructor: Dr. Khalid Raza (JMI)
# Assignment Q1–Q6 coverage per brief (dataset → mapping/QC → quant/normalization
# → DEGs → functional enrichment → optional ML gene ranking)
# ----------------------------------------------------------------------------
# Notes
# - Uses TCGA-BRCA RNA-seq (STAR counts) via TCGAbiolinks.
# - Saves everything under ./TCGA_BRCA_Analysis
# - Robust to list-like metadata, missing enrichment hits, and SVM prob preds.
# - Prefers TPM provided by TCGA (tpm_unstrand) over ad‑hoc FPKM lengths.
# - All figures are 300 DPI PNGs; a single HTML report is generated.
# ============================================================================

# ----------------------------- 0) Environment -------------------------------
options(stringsAsFactors = FALSE, scipen = 999, max.print = 200)
set.seed(1234)
gc()

cat("\n========== INITIALIZING CLEAN PIPELINE ==========%n")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Required packages (Bioc + CRAN)
required_packages <- c(
  # Bioconductor
  "TCGAbiolinks", "DESeq2", "edgeR", "limma", "SummarizedExperiment",
  "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler", "enrichplot", "DOSE",
  # CRAN
  "tidyverse", "data.table", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
  "cowplot", "vsn", "VennDiagram", "caret", "randomForest", "glmnet", "pROC",
  "viridis"
)

cat("Checking/Installing packages...\n")
for (pkg in required_packages) {
  suppressPackageStartupMessages({
    if (!suppressWarnings(require(pkg, character.only = TRUE))) {
      if (pkg %in% c("tidyverse","data.table","ggplot2","ggrepel","pheatmap",
                     "RColorBrewer","cowplot","VennDiagram","caret",
                     "randomForest","glmnet","pROC","viridis")) {
        install.packages(pkg, dependencies = TRUE)
      } else {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      }
      library(pkg, character.only = TRUE)
    }
  })
}

# Dir scaffold
base_dir <- "TCGA_BRCA_Analysis"
dirs <- file.path(base_dir, c("data","results","figures","reports","QC","DEGs","enrichment","ML"))
dir.create(base_dir, showWarnings = FALSE)
for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Helpers --------------------------------------------------------------------
flatten_df <- function(df) {
  df <- as.data.frame(df)
  is_listcol <- vapply(df, function(x) is.list(x) || inherits(x, "List") || inherits(x, "CompressedList"), logical(1))
  if (any(is_listcol)) {
    df[is_listcol] <- lapply(df[is_listcol], function(col) {
      vapply(col, function(y) paste(na.omit(as.character(y)), collapse = "; "), character(1L))
    })
  }
  df
}

safe_write_csv <- function(x, path, row.names = FALSE) {
  tryCatch({ utils::write.csv(x, path, row.names = row.names) }, error = function(e) {
    utils::write.csv(flatten_df(x), path, row.names = row.names)
  })
}

open_png <- function(path, w = 12, h = 8, dpi = 300, units = "in") {
  png(path, width = w, height = h, res = dpi, units = units)
}

# ----------------------------- Q1) Dataset ----------------------------------
cat("\n========== Q1: DATASET ==========%n")
cat("Project: TCGA-BRCA | Data: STAR Counts + TPM (if available)\n")

query_tumor <- TCGAbiolinks::GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  experimental.strategy = "RNA-Seq"
)

query_normal <- TCGAbiolinks::GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal",
  experimental.strategy = "RNA-Seq"
)

tumor_meta  <- TCGAbiolinks::getResults(query_tumor)
normal_meta <- TCGAbiolinks::getResults(query_normal)
cat(sprintf("Tumor available: %d | Normal available: %d\n", nrow(tumor_meta), nrow(normal_meta)))

# Memory-friendly selection (≤100 tumor + all normals)
set.seed(123)
max_tumor <- min(100, nrow(tumor_meta))
sel_tumor  <- sample(tumor_meta$cases, max_tumor)
sel_normal <- normal_meta$cases

query_final <- TCGAbiolinks::GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(sel_tumor, sel_normal)
)

cat("Downloading (GDC API)...\n")
tryCatch({
  TCGAbiolinks::GDCdownload(query_final, method = "api", files.per.chunk = 10)
}, error = function(e) {
  cat("API failed, trying GDC client...\n"); TCGAbiolinks::GDCdownload(query_final, method = "client")
})

cat("Preparing SummarizedExperiment...\n")
se <- TCGAbiolinks::GDCprepare(query_final, summarizedExperiment = TRUE)

# Save raw object
saveRDS(se, file.path(base_dir, "data/tcga_brca_raw_data.rds"))
cat("✓ Data prepared and saved.\n")

# Sample metadata (flattened)
sample_info <- flatten_df(SummarizedExperiment::colData(se))
safe_write_csv(sample_info, file.path(base_dir, "data/sample_metadata.csv"))

# Summary
cat(sprintf("Samples total: %d (Tumor=%d, Normal=%d)\n",
            ncol(se), sum(sample_info$sample_type == "Primary Tumor"), sum(sample_info$sample_type == "Solid Tissue Normal")))
cat("Assays detected: ", paste(SummarizedExperiment::assayNames(se), collapse=", "), "\n")

# ----------------------------- Q2) Mapping/QC -------------------------------
cat("\n========== Q2: READ MAPPING & QC ==========%n")
cat("TCGA pipeline uses STAR (two-pass) on GRCh38 with GENCODE annotation.\n")

# Basic QC visuals using counts assay
counts_assay_name <- if ("unstranded" %in% assayNames(se)) "unstranded" else assayNames(se)[1]
counts <- SummarizedExperiment::assay(se, counts_assay_name)

open_png(file.path(base_dir, "figures/Q2_alignment_summary.png"), w=12, h=10)
par(mfrow=c(2,2))
barplot(table(sample_info$sample_type), col=c("#E74C3C","#3498DB"), main="Sample Distribution", ylab="Count", las=2)
# Genes detected per sample
gene_detect <- colSums(counts > 0)
hist(gene_detect, breaks = 50, col = "skyblue", main = "Genes Detected per Sample", xlab = "+ counts > 0")
# Library size per group
lib_size <- colSums(counts)
boxplot(lib_size ~ sample_info$sample_type, col=c("#3498DB","#E74C3C"), main="Library Size", ylab="Total counts", las=2)
# Mean expression density
mean_expr <- rowMeans(counts)
plot(density(log2(mean_expr + 1)), main="Mean Expression (log2)", xlab="log2(mean+1)", col="darkblue", lwd=2)
dev.off()

safe_write_csv(data.frame(barcode = colnames(se), lib_size), file.path(base_dir, "QC/library_sizes.csv"))

# ---------------------- Q3) Quantification & Normalization ------------------
cat("\n========== Q3: EXPRESSION QUANTIFICATION & NORMALIZATION ==========%n")

# Filter low expression (counts >10 in ≥10% samples)
min_samples <- ceiling(0.10 * ncol(counts))
keep <- rowSums(counts > 10) >= min_samples
counts_f <- counts[keep, , drop=FALSE]
cat(sprintf("Genes (raw)=%d | kept=%d | removed=%d\n", nrow(counts), nrow(counts_f), nrow(counts)-nrow(counts_f)))

# Sample metadata aligned to columns
condition <- ifelse(sample_info$sample_type == "Primary Tumor", "Tumor", "Normal")
coldata <- data.frame(condition = factor(condition, levels = c("Normal","Tumor")), row.names = colnames(counts_f))

# DESeq2 object + size factors
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_f, colData = coldata, design = ~ condition)
dds <- DESeq2::estimateSizeFactors(dds)
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
safe_write_csv(norm_counts, file.path(base_dir, "data/normalized_counts.csv"))

# Prefer TPM if present
TPM <- NULL
if ("tpm_unstrand" %in% assayNames(se)) {
  TPM <- SummarizedExperiment::assay(se, "tpm_unstrand")[rownames(counts_f), colnames(counts_f), drop=FALSE]
  safe_write_csv(TPM, file.path(base_dir, "data/tpm_counts.csv"))
}

# QC plots (boxplots, size factors, PCA)
open_png(file.path(base_dir, "figures/Q3_normalization_QC.png"), w=14, h=10)
par(mfrow=c(2,3))
boxplot(log2(counts[ , 1:min(20,ncol(counts))] + 1), col="lightblue", main="Raw Counts (first 20)", ylab="log2(Counts+1)", las=2, cex.axis=0.7)
boxplot(log2(norm_counts[, 1:min(20,ncol(norm_counts))] + 1), col="lightgreen", main="Normalized Counts", ylab="log2(Norm+1)", las=2, cex.axis=0.7)
plot(DESeq2::sizeFactors(dds), pch=19, col=ifelse(coldata$condition=="Tumor","red","blue"), main="DESeq2 Size Factors", ylab="Size factor"); legend("topright", c("Tumor","Normal"), col=c("red","blue"), pch=19)
vsd <- DESeq2::vst(dds, blind = FALSE)
print( DESeq2::plotPCA(vsd, intgroup="condition") + ggplot2::ggtitle("PCA – VST") + ggplot2::theme_minimal() )
# Mean-SD plot (vsn)
vsn::meanSdPlot(SummarizedExperiment::assay(vsd))
plot(stats::density(colSums(norm_counts)), main="Density of Library Sizes (normalized)", xlab="Sum of normalized counts")
dev.off()

norm_summary <- data.frame(
  Metric = c("Samples","Tumor","Normal","Genes(raw)","Genes(kept)","Median lib Tumor","Median lib Normal","Mean size factor"),
  Value  = c(ncol(counts_f), sum(coldata$condition=="Tumor"), sum(coldata$condition=="Normal"), nrow(counts), nrow(counts_f),
             median(colSums(counts_f[, coldata$condition=="Tumor", drop=FALSE])),
             median(colSums(counts_f[, coldata$condition=="Normal", drop=FALSE])),
             mean(DESeq2::sizeFactors(dds)))
)
safe_write_csv(norm_summary, file.path(base_dir, "results/Q3_normalization_summary.csv"))

# --------------------------- Q4) DEG Analysis -------------------------------
cat("\n========== Q4: DEGs (DESeq2) ==========%n")

dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c("condition","Tumor","Normal"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Add gene symbols
ens <- sub("\\..*","", res_df$gene_id)
ann <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ens), keytype = "ENSEMBL", columns = c("ENSEMBL","SYMBOL","GENENAME"))
res_df$gene_symbol <- ann$SYMBOL[match(ens, ann$ENSEMBL)]

safe_write_csv(res_df, file.path(base_dir, "DEGs/all_deseq2_results.csv"))

# Assignment thresholds
DEG <- res_df %>% dplyr::filter(abs(log2FoldChange) >= 2 & pvalue <= 0.05) %>% dplyr::arrange(pvalue)
safe_write_csv(DEG, file.path(base_dir, "DEGs/significant_DEGs.csv"))

cat(sprintf("Genes tested=%d | Significant (|log2FC|>=2, p<=0.05)=%d\n", nrow(res_df), nrow(DEG)))

# Build plotting df with status flag
plot_df <- res_df %>%
  dplyr::mutate(
    status = dplyr::case_when(
      log2FoldChange >=  2 & pvalue <= 0.05 ~ "Up",
      log2FoldChange <= -2 & pvalue <= 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Choose up to 15 labelled points from significant (non-NS), non-NA symbols
label_df <- plot_df %>%
  dplyr::filter(status != "NS", !is.na(gene_symbol)) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::slice_head(n = 15)

print(
  ggplot(plot_df, aes(log2FoldChange, -log10(pvalue), color = status)) +
    geom_point(alpha = .6, size = 1.2) +
    scale_color_manual(values = c(Up = "#E74C3C", Down = "#3498DB", NS = "grey70")) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene_symbol),
      inherit.aes = TRUE,
      size = 3,
      max.overlaps = 15
    ) +
    theme_minimal() +
    labs(title = "Volcano: Tumor vs Normal", x = "log2FC", y = "-log10(p)")
)

# Boxplots top5
top5 <- head(DEG, 5)
open_png(file.path(base_dir, "figures/Q4_top5_DEGs_boxplots.png"), w=14, h=10)
par(mfrow=c(2,3))
for (i in seq_len(nrow(top5))) {
  g <- top5$gene_id[i]; sym <- top5$gene_symbol[i]
  y <- norm_counts[g, ]
  boxplot(y ~ coldata$condition, col=c("#3498DB","#E74C3C"), main=paste0(sym, " (", g, ")"), ylab="Normalized expression", las=2)
  stripchart(y ~ coldata$condition, vertical=TRUE, method="jitter", pch=20, add=TRUE)
}
dev.off()

# MA plot
open_png(file.path(base_dir, "figures/Q4_MA_plot.png"), w=10, h=8)
DESeq2::plotMA(res, ylim=c(-5,5), main="MA Plot")
dev.off()

# Heatmap top50
top50_ids <- head(DEG$gene_id, 50)
mat50 <- norm_counts[top50_ids, , drop=FALSE]
rownames(mat50) <- DEG$gene_symbol[match(rownames(mat50), DEG$gene_id)]
mat50s <- t(scale(t(mat50)))
ann_col <- data.frame(Condition = coldata$condition, row.names = colnames(mat50s))
ann_colors <- list(Condition = c(Tumor="#E74C3C", Normal="#3498DB"))
open_png(file.path(base_dir, "figures/Q4_heatmap_top50_DEGs.png"), w=12, h=14)
pheatmap::pheatmap(mat50s, annotation_col = ann_col, annotation_colors = ann_colors,
                   color = colorRampPalette(c("blue","white","red"))(100), show_colnames=FALSE,
                   fontsize_row=8, main="Top 50 DEGs (scaled)")
dev.off()

# ---------------------- Q5) Functional Enrichment ---------------------------
cat("\n========== Q5: FUNCTIONAL ENRICHMENT ==========%n")

# Gene symbol → Entrez
deg_symbols <- DEG$gene_symbol[!is.na(DEG$gene_symbol)]
if (length(deg_symbols) == 0) stop("No gene symbols available for enrichment.")

conv <- clusterProfiler::bitr(deg_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ents <- unique(conv$ENTREZID)

# GO BP/MF/CC
GO_list <- list()
for (ont in c("BP","MF","CC")) {
  eg <- tryCatch(enrichGO(gene = ents, OrgDb = org.Hs.eg.db, ont = ont, pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE), error=function(e) NULL)
  GO_list[[ont]] <- eg
  if (!is.null(eg) && nrow(as.data.frame(eg)) > 0) {
    safe_write_csv(as.data.frame(eg), file.path(base_dir, sprintf("enrichment/GO_%s.csv", ont)))
  }
}

open_png(file.path(base_dir, "figures/Q5_GO_enrichment.png"), w=14, h=12)
plots <- list()
for (ont in c("BP","MF","CC")) {
  obj <- GO_list[[ont]]
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    plots[[ont]] <- enrichplot::dotplot(obj, showCategory = 15, title = paste("GO",ont))
  } else {
    plots[[ont]] <- ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::labs(title=paste("GO",ont, "(no hits)"))
  }
}
print(cowplot::plot_grid(plots$BP, plots$MF, plots$CC, ncol = 1))
dev.off()

# KEGG
kegg <- tryCatch(enrichKEGG(gene = ents, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2), error=function(e) NULL)
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  safe_write_csv(as.data.frame(kegg), file.path(base_dir, "enrichment/KEGG_pathways.csv"))
  open_png(file.path(base_dir, "figures/Q5_KEGG_pathways.png"), w=12, h=10)
  print(enrichplot::dotplot(kegg, showCategory = 20, title = "KEGG Pathways"))
  dev.off()
} else {
  cat("No significant KEGG pathways.\n")
}

# Disease Ontology
edo <- tryCatch(enrichDO(gene = ents, ont = "DO", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE), error=function(e) NULL)
if (!is.null(edo) && nrow(as.data.frame(edo)) > 0) {
  safe_write_csv(as.data.frame(edo), file.path(base_dir, "enrichment/DiseaseOntology.csv"))
  open_png(file.path(base_dir, "figures/Q5_disease_ontology.png"), w=12, h=10)
  print(enrichplot::barplot(edo, showCategory = 20, title = "Disease Ontology"))
  dev.off()
}

# NCG (optional; may not be available in all setups)
ncg <- tryCatch(enrichNCG(gene = ents, pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE), error=function(e) NULL)
if (!is.null(ncg) && nrow(as.data.frame(ncg)) > 0) {
  safe_write_csv(as.data.frame(ncg), file.path(base_dir, "enrichment/NCG_cancer_genes.csv"))
  open_png(file.path(base_dir, "figures/Q5_cancer_genes.png"), w=12, h=10)
  print(enrichplot::dotplot(ncg, showCategory = 20, title = "Network of Cancer Genes"))
  dev.off()
}

# ---------------- Q6) Machine Learning – Gene Ranking (Optional) ------------
cat("\n========== Q6: MACHINE LEARNING (OPTIONAL) ==========%n")

if (nrow(DEG) >= 10) {
  top10 <- head(DEG, 10)
  top10_ids <- top10$gene_id
  gene_names <- top10$gene_symbol[match(top10_ids, top10$gene_id)]
  
  # Feature matrix: log2(TPM+1) preferred; fallback to log2(norm+1)
  if (!is.null(TPM)) {
    X <- log2(TPM[top10_ids, , drop=FALSE] + 1)
  } else {
    X <- log2(norm_counts[top10_ids, , drop=FALSE] + 1)
  }
  X <- t(X)
  colnames(X) <- gene_names
  y <- coldata$condition
  
  ml_df <- data.frame(X, condition = y)
  safe_write_csv(ml_df, file.path(base_dir, "ML/top10_features.csv"), row.names = TRUE)
  
  # Train/test split
  set.seed(42)
  idx <- caret::createDataPartition(ml_df$condition, p = 0.7, list = FALSE)
  tr <- ml_df[idx, ]; te <- ml_df[-idx, ]
  
  # Random Forest
  rf <- randomForest::randomForest(condition ~ ., data = tr, ntree = 500, importance = TRUE, mtry = max(1, floor(sqrt(ncol(tr)-1))))
  rf_pred <- predict(rf, te)
  rf_conf <- caret::confusionMatrix(rf_pred, te$condition)
  rf_imp <- as.data.frame(randomForest::importance(rf))
  rf_imp$Gene <- rownames(rf_imp)
  rf_imp <- rf_imp[order(rf_imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
  safe_write_csv(rf_imp, file.path(base_dir, "ML/rf_feature_importance.csv"))
  
  open_png(file.path(base_dir, "figures/Q6_RF_feature_importance.png"), w=10, h=8)
  randomForest::varImpPlot(rf, main = "Random Forest – Feature Importance")
  dev.off()
  
  # LASSO Logistic
  x_tr <- as.matrix(tr[, setdiff(names(tr), "condition")]); y_tr <- tr$condition
  x_te <- as.matrix(te[, setdiff(names(te), "condition")]); y_te <- te$condition
  cv <- glmnet::cv.glmnet(x_tr, y_tr, family="binomial", alpha=1, type.measure = "class")
  lasso <- glmnet::glmnet(x_tr, y_tr, family="binomial", alpha=1, lambda = cv$lambda.min)
  lprob <- as.numeric(predict(lasso, x_te, type = "response"))
  lpred <- factor(ifelse(lprob > 0.5, "Tumor","Normal"), levels=c("Normal","Tumor"))
  l_conf <- caret::confusionMatrix(lpred, y_te)
  lcoef <- as.matrix(coef(lasso))
  ldf <- data.frame(Gene = rownames(lcoef), Coefficient = lcoef[,1]) %>% dplyr::filter(Gene!="(Intercept)") %>% dplyr::arrange(dplyr::desc(abs(Coefficient)))
  safe_write_csv(ldf, file.path(base_dir, "ML/lasso_coefficients.csv"))
  
  open_png(file.path(base_dir, "figures/Q6_LASSO_coefficients.png"), w=10, h=8)
  print( ggplot(ldf, aes(x = reorder(Gene, abs(Coefficient)), y = Coefficient)) + geom_bar(stat="identity", fill="#3498DB", alpha=.85) + coord_flip() + theme_minimal() + labs(title="LASSO Coefficients", x="Gene") )
  dev.off()
  
  # SVM (with probabilities)
  ctrl <- caret::trainControl(method = "cv", number = 5, classProbs = TRUE)
  svmM <- caret::train(condition ~ ., data = tr, method = "svmRadial", trControl = ctrl, preProcess = c("center","scale"))
  spred <- predict(svmM, te)
  s_conf <- caret::confusionMatrix(spred, te$condition)
  
  # ROC curves
  open_png(file.path(base_dir, "figures/Q6_ROC_curves.png"), w=12, h=8)
  rf_prob <- predict(rf, te, type = "prob")[, "Tumor"]
  s_prob  <- predict(svmM, te, type = "prob")[, "Tumor"]
  rroc <- pROC::roc(te$condition, rf_prob, levels=c("Normal","Tumor"), direction="<")
  lroc <- pROC::roc(te$condition, lprob,  levels=c("Normal","Tumor"), direction="<")
  sroc <- pROC::roc(te$condition, s_prob, levels=c("Normal","Tumor"), direction="<")
  plot(rroc, col="#E74C3C", lwd=2, main="ROC – Model Comparison")
  plot(lroc, col="#3498DB", lwd=2, add=TRUE)
  plot(sroc, col="#2ECC71", lwd=2, add=TRUE)
  legend("bottomright", legend = c(paste0("RF (AUC=", round(pROC::auc(rroc),3), ")"), paste0("LASSO (AUC=", round(pROC::auc(lroc),3), ")"), paste0("SVM (AUC=", round(pROC::auc(sroc),3), ")")), col=c("#E74C3C","#3498DB","#2ECC71"), lwd=2)
  dev.off()
  
  # Comparison table
  comp <- data.frame(
    Model = c("Random Forest","LASSO","SVM"),
    Accuracy    = c(rf_conf$overall["Accuracy"], l_conf$overall["Accuracy"], s_conf$overall["Accuracy"]),
    Sensitivity = c(rf_conf$byClass["Sensitivity"], l_conf$byClass["Sensitivity"], s_conf$byClass["Sensitivity"]),
    Specificity = c(rf_conf$byClass["Specificity"], l_conf$byClass["Specificity"], s_conf$byClass["Specificity"]),
    AUC         = c(as.numeric(pROC::auc(rroc)), as.numeric(pROC::auc(lroc)), as.numeric(pROC::auc(sroc)))
  )
  safe_write_csv(comp, file.path(base_dir, "ML/model_comparison.csv"))
  
  # Combined gene ranking (RF + LASSO)
  rank <- data.frame(Gene = rf_imp$Gene,
                     RF_Importance = rf_imp$MeanDecreaseAccuracy,
                     LASSO_Coefficient = abs(ldf$Coefficient[match(rf_imp$Gene, ldf$Gene)]))
  rank$Combined_Score <- as.numeric(scale(rank$RF_Importance)) + as.numeric(scale(rank$LASSO_Coefficient))
  rank <- rank[order(rank$Combined_Score, decreasing = TRUE), ]
  safe_write_csv(rank, file.path(base_dir, "ML/gene_ranking_combined.csv"))
  
  # Heatmap of importances
  open_png(file.path(base_dir, "figures/Q6_gene_ranking_heatmap.png"), w=10, h=8)
  m <- as.matrix(rank[, c("RF_Importance","LASSO_Coefficient")]); rownames(m) <- rank$Gene
  pheatmap::pheatmap(scale(m), cluster_rows = TRUE, cluster_cols = FALSE, color = colorRampPalette(c("white","orange","red"))(100), main="Gene Importance (RF + LASSO)")
  dev.off()
  
  cat("✓ ML completed.\n")
} else {
  cat("<10 DEGs – skipping ML section.\n")
}

# ---------------------- Report & Session Info -------------------------------
cat("\n========== REPORTS ==========%n")

# Simple HTML report (concise)
html <- file.path(base_dir, "reports/TCGA_BRCA_Report.html")
cat("Generating HTML report...\n")

report <- list()
report$header <- sprintf("<h1>TCGA-BRCA RNA-seq Analysis</h1><p>Date: %s</p>", format(Sys.time(), "%Y-%m-%d %H:%M"))
report$q1 <- sprintf("<h2>Q1 Dataset</h2><ul><li>Samples: %d (Tumor=%d, Normal=%d)</li><li>Assays: %s</li><li>Counts assay: %s</li></ul>", ncol(se), sum(coldata$condition=="Tumor"), sum(coldata$condition=="Normal"), paste(assayNames(se), collapse=", "), counts_assay_name)
report$q2 <- "<h2>Q2 Mapping & QC</h2><p>STAR aligner on GRCh38 (two-pass). See QC figures.</p><img src='../figures/Q2_alignment_summary.png' width='900'>"
report$q3 <- sprintf("<h2>Q3 Normalization</h2><ul><li>Filter: counts>10 in ≥10%% samples</li><li>Method: DESeq2 size factors (%0.3f mean)</li><li>TPM preferred: %s</li></ul><img src='../figures/Q3_normalization_QC.png' width='900'>",
                     mean(DESeq2::sizeFactors(dds)), ifelse(!is.null(TPM),"yes","no"))
report$q4 <- sprintf("<h2>Q4 DEGs</h2><ul><li>Tested: %d</li><li>Significant (|log2FC|≥2, p≤0.05): <b>%d</b> (Up=%d, Down=%d)</li></ul><img src='../figures/Q4_volcano_plot.png' width='900'><img src='../figures/Q4_top5_DEGs_boxplots.png' width='900'><img src='../figures/Q4_heatmap_top50_DEGs.png' width='900'>",
                     nrow(res_df), nrow(DEG), sum(DEG$log2FoldChange>2), sum(DEG$log2FoldChange < -2))
report$q5 <- "<h2>Q5 Enrichment</h2><p>See GO/KEGG/DO figures and CSVs under enrichment/.</p><img src='../figures/Q5_GO_enrichment.png' width='900'>"
if (file.exists(file.path(base_dir, "figures/Q5_KEGG_pathways.png"))) report$q5 <- paste0(report$q5, "<img src='../figures/Q5_KEGG_pathways.png' width='900'>")
if (file.exists(file.path(base_dir, "figures/Q5_disease_ontology.png"))) report$q5 <- paste0(report$q5, "<img src='../figures/Q5_disease_ontology.png' width='900'>")

if (file.exists(file.path(base_dir, "ML/model_comparison.csv"))) {
  comp <- read.csv(file.path(base_dir, "ML/model_comparison.csv"))
  report$q6 <- sprintf("<h2>Q6 ML (Optional)</h2><p>Models trained: RF/LASSO/SVM.</p><pre>%s</pre><img src='../figures/Q6_ROC_curves.png' width='900'><img src='../figures/Q6_RF_feature_importance.png' width='700'><img src='../figures/Q6_LASSO_coefficients.png' width='700'><img src='../figures/Q6_gene_ranking_heatmap.png' width='700'>", paste(capture.output(print(comp, row.names=FALSE)), collapse = "\n"))
} else {
  report$q6 <- "<h2>Q6 ML (Optional)</h2><p>Skipped (insufficient DEGs).</p>"
}

html_txt <- paste0("<html><head><meta charset='utf-8'><style>body{font-family:Arial;margin:30px} img{border:1px solid #ccc;margin:10px 0} pre{background:#f6f8fa;padding:10px;border:1px solid #e1e4e8;overflow:auto}</style></head><body>",
                   report$header, report$q1, report$q2, report$q3, report$q4, report$q5, report$q6,
                   "<hr><p><i>Generated by Clean Pipeline (R).</i></p></body></html>")
writeLines(html_txt, html)
cat("✓ HTML report:", html, "\n")

# Text summary
sink(file.path(base_dir, "reports/Analysis_Summary.txt"))
cat("TCGA-BRCA RNA-seq Analysis – Summary\n===============================\n\n")
print(norm_summary)
cat("\nDEG counts:\n"); cat(sprintf("Tested=%d | Significant=%d | Up=%d | Down=%d\n", nrow(res_df), nrow(DEG), sum(DEG$log2FoldChange>2), sum(DEG$log2FoldChange< -2)))
if (exists("comp")) { cat("\nModel comparison:\n"); print(comp) }
sink()

# Session info
writeLines(capture.output(sessionInfo()), file.path(base_dir, "reports/session_info.txt"))
cat("✓ Session info saved.\n")

cat("\n████ CLEAN PIPELINE FINISHED SUCCESSFULLY ████\n")
cat(sprintf("Output folder: %s\n", base_dir))
cat("Key files:\n - data/normalized_counts.csv\n - DEGs/significant_DEGs.csv\n - figures/Q4_volcano_plot.png\n - enrichment/GO_*.csv, KEGG_pathways.csv, DiseaseOntology.csv\n - reports/TCGA_BRCA_Report.html\n")
