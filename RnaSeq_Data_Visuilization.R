# RNA-seq Differential Expression Analysis Script
# Author: Parth Doshi
# Description: This script analyzes differential expression of genes from GEO dataset GSE37764 and visualizes results.
# Dependencies: DESeq2, AnnotationDbi, org.Hs.eg.db, GEOquery, ggplot2, pheatmap, EnhancedVolcano, WGCNA

# 1. Load and Install Required Packages ---------------------------------------------------------

required_packages <- c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "GEOquery", 
                       "ggplot2", "pheatmap", "EnhancedVolcano", "WGCNA")

install_and_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message(paste("Installing package:", package))
    if (package %in% c("DESeq2", "org.Hs.eg.db", "GEOquery", "WGCNA")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(package)
    } else {
      install.packages(package)
    }
  }
  library(package, character.only = TRUE)
}

invisible(lapply(required_packages, install_and_load))

# Enable multi-threading in WGCNA
enableWGCNAThreads()

# 2. Define Parameters ---------------------------------------------------------------------------

geo_id <- "GSE37764"
counts_file <- "C:/Users/Parth Doshi/Desktop/Rna-seq_DE_gene_CO_expression_graph_theory/GSE37764_raw_counts_GRCh38.p13_NCBI.tsv.gz"

# 3. Load GEO Metadata and Prepare Clinical Data -------------------------------------------------

geo_data <- getGEO(geo_id, GSEMatrix = TRUE)
clin_data <- pData(phenoData(geo_data[[1]]))

if ("tissue type:ch1" %in% colnames(clin_data)) {
  clin_data$condition <- factor(clin_data$`tissue type:ch1`)
} else if ("characteristics_ch1" %in% colnames(clin_data)) {
  clin_data$condition <- factor(clin_data$characteristics_ch1)
} else {
  stop("Cannot find a column for 'condition' in the GEO metadata.")
}
message(paste("Loaded metadata for", nrow(clin_data), "samples."))

# 4. Load and Subset Raw Count Data --------------------------------------------------------------

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
common_samples <- intersect(colnames(counts), rownames(clin_data))
if (length(common_samples) == 0) stop("No matching samples found between count data and metadata.")
counts <- counts[, common_samples]
clin_data <- clin_data[common_samples, ]

# 5. Differential Expression Analysis with DESeq2 ------------------------------------------------

coldata <- data.frame(condition = clin_data$condition)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumor", "normal"))
res <- res[order(res$padj), ]

# 6. Visualization of Differential Expression Results --------------------------------------------

# Volcano Plot
res$GeneSymbol <- mapIds(org.Hs.eg.db, keys = rownames(res), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
volcano <- EnhancedVolcano(res,
                           lab = res$GeneSymbol,
                           x = 'log2FoldChange',
                           y = 'padj',
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           pointSize = 3.0,
                           labSize = 4.0,
                           title = 'Differential Expression: Tumor vs Normal',
                           subtitle = 'Volcano Plot')
ggsave(filename = "volcano_plot_tumor_vs_normal.png", plot = volcano, width = 8, height = 6)

# PCA Plot
rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
ggsave(filename = "pca_plot.png", plot = pca_plot, width = 8, height = 6)

# Heatmap of Top Differentially Expressed Genes
vsd <- vst(dds, blind = FALSE)
topGenes <- head(order(res$padj), 50)
heatmap_data <- assay(vsd)[topGenes, ]
heatmap_data <- heatmap_data - rowMeans(heatmap_data)
heatmap_plot <- pheatmap(heatmap_data,
                         annotation_col = as.data.frame(colData(vsd)[, "condition", drop = FALSE]),
                         show_rownames = FALSE,
                         main = "Top 50 Differentially Expressed Genes")
ggsave(filename = "heatmap.png", plot = heatmap_plot, width = 10, height = 8)

# MA Plot
ma_plot <- plotMA(res, ylim = c(-5, 5), main = "MA Plot")
ggsave(filename = "ma_plot.png", plot = ma_plot, width = 8, height = 6)

# 7. WGCNA Analysis ------------------------------------------------------------------------------

# Filter for High Variance Genes
variance <- apply(assay(vsd), 1, var)
expr_data <- assay(vsd)[variance > quantile(variance, 0.75), ]
rownames(expr_data) <- rownames(assay(vsd))[variance > quantile(variance, 0.75)]
expr_data <- t(expr_data)

# Soft Thresholding Power Selection
powers <- c(1:20)
sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 5)
softPower <- 8

# Adjacency and TOM Calculation
adjacency <- adjacency(expr_data, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Hierarchical Clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamic Tree Cut for Module Detection
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)

# Dendrogram Plot with Module Colors
pdf("dendrogram_all_high_variance_genes.pdf", width = 12, height = 8)
plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Dendrogram with Dynamic Tree Cut Colors (All High-Variance Genes)"
)
dev.off()

# 8. Module-Trait Relationship Analysis ----------------------------------------------------------

design <- model.matrix(~ condition - 1, data = clin_data)
colnames(design) <- levels(clin_data$condition)

MEs <- moduleEigengenes(expr_data, colors = dynamicColors)$eigengenes
moduleTraitCor <- cor(MEs, design, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(expr_data))

# Heatmap of Module-Trait Relationships
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

png("module_trait_relationships_conditions_no_labels_adjusted.png", width = 1500, height = 1800, res = 150)
par(mar = c(6, 6, 4, 2))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(design),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = NULL,
  setStdMargins = TRUE,
  cex.lab = 0.8,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships: Tumor vs Normal",
  xLabelsAngle = 45,
  xLabelsAdj = 1
)

dev.off()
