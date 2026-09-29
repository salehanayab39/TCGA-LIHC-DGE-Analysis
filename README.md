 TCGA-LIHC Differential Gene Expression Analysis

Overview
This project performs a differential gene expression (DGE) analysis of 
liver cancer (hepatocellular carcinoma) versus normal liver tissue using 
public RNA-sequencing data from The Cancer Genome Atlas (TCGA-LIHC).

 Research Question
Which genes are significantly differentially expressed between tumour and 
normal liver tissue, and what biological processes do these genes regulate?

 Data
-Source: TCGA-LIHC via Genomic Data Commons (GDC)
-Samples:374 primary tumours + 50 solid tissue normals (424 total)
-Data type: RNA-seq, STAR-Counts workflow
-Genes analysed:49,319 (after filtering low-count genes)

 Methods
- Data download: TCGAbiolinks` R package
- Differential expression: `DESeq2` (padj < 0.05, |log2FC| > 1)
- Visualisation: `EnhancedVolcano`, `pheatmap`
- Functional enrichment:`clusterProfiler` (GO Biological Process)

 Key Results
- 9,457 genes significantly differentially expressed
- 7,529 upregulated in tumour tissue
- 1,928 downregulated in tumour tissue

 Top Upregulated Genes

Gene	Log2FC	padj	Function		
AFP	+6.00	2.2e-39	Fetal marker / diagnostic		
GPC3	+5.92	3.5e-80	HCC cell surface marker		
ZIC2	+6.70	4.9e-86	Stemness transcription factor		
CDKN3	+3.90	1.6e-90	Cell cycle phosphatase		
BIRC5	+4.10	1.3e-82	Anti-apoptosis (Survivin)		


 Known HCC Markers Rediscovered
- AFP (up +6.00) — clinical blood biomarker
- GPC3 (up +5.92) — immunohistochemistry marker
- ALB (down −1.84) — liver function lost
- APOA1 (down −1.33) — lipid metabolism lost
- SPP1 (up +4.30) — inflammation/fibrosis

 Biological Theme
Upregulated genes were dominated by cell cycle, mitosis, and 
DNA replication pathways — consistent with the proliferative 
nature of hepatocellular carcinoma.
                                                                                                                                                                   Repository Contents                                                                                                                                                    ├── analysis.R # Complete analysis script with comments
├── figures/
│ ├── volcano_plot.png # Volcano plot of all genes
│ ├── heatmap.png # Top 50 genes heatmap
│ └── go_enrichment.png # GO enrichment dot plot
├── report/
│ └── TCGA_LIHC_Report.pdf # 2-3 page written report
└── results/
└── TCGA_LIHC_DESeq2_results.csv # Full DESeq2 results
                                                                                                                                                                               Tools Used
- R (version 4.6.1)
- RStudio
- Packages: TCGAbiolinks, DESeq2, EnhancedVolcano, pheatmap, clusterProfiler, org.Hs.eg.db

 Author
Saleha Nayab
salehanayab39@gmail.com

 License
This project is open source and available under the MIT License.