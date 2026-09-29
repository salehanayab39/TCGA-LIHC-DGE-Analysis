# ============================================
# TCGA-LIHC Differential Gene Expression Analysis
# Author: Saleha Nayab
# 

# ---- SETUP ----
rm(list = ls())
options(timeout = 2400)
setwd("C:/Users/saleh/Desktop/TCGA_LIHC_Project")

# ---- LOAD LIBRARIES ----
library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# ---- QUERY TCGA-LIHC DATA ----
query <- GDCquery(
  project = 'TCGA-LIHC',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

# ---- DOWNLOAD DATA (run once) ----
# GDCdownload(query, method = "api", files.per.chunk = 5)

# ---- PREPARE DATA ----
data <- GDCprepare(query)

# ---- DIFFERENTIAL EXPRESSION ANALYSIS ----
counts <- assay(data, 'unstranded')
coldata <- colData(data)
coldata$condition <- ifelse(grepl('Tumor', coldata$sample_type), 'tumour', 'normal')
keep <- !is.na(coldata$condition)
counts <- counts[, keep]
coldata <- coldata[keep, ]

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c('condition', 'tumour', 'normal'))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Add gene symbols
res_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = sub("\\..*", "", res_df$gene),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# ---- VOLCANO PLOT ----
EnhancedVolcano(res_df,
                lab = res_df$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Liver Cancer vs Normal',
                subtitle = 'TCGA-LIHC'
)

# ---- HEATMAP ----
sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
sig_genes <- sig_genes[order(sig_genes$padj), ]
top50 <- head(sig_genes, 50)

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top50$gene, ]
mat_scaled <- t(scale(t(mat)))
annotation <- data.frame(condition = coldata$condition, row.names = colnames(mat))

pheatmap(mat_scaled,
         annotation_col = annotation,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = 'Top 50 Differentially Expressed Genes'
)

# ---- GO ENRICHMENT ----
up_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1, ]
gene_ids <- bitr(sub("\\..*", "", up_genes$gene),
                 fromType = 'ENSEMBL',
                 toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db
)

ego <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego, showCategory = 20, title = 'GO Biological Process Enrichment')

# ---- SAVE RESULTS ----
write.csv(res_df, "results/TCGA_LIHC_DESeq2_results.csv", row.names = FALSE)