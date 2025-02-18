### RNA-seq Differential Expression and Co-expression Analysis Pipeline

#### Project Overview:
This project is a comprehensive RNA-seq analysis pipeline designed to perform **differential expression analysis**, **co-expression network analysis** (using WGCNA), and various **data visualizations**. The analysis is based on publicly available **GEO datasets** and uses RNA-seq raw count data. It integrates key steps in the analysis of transcriptomic data, such as normalization, differential expression testing, and co-expression analysis, with clear visual outputs for biological interpretation.

---

#### Key Features:
- **Automated RNA-seq Pipeline**: The entire process from loading GEO datasets to generating publication-quality visualizations is fully automated.
- **Differential Expression Analysis**: The pipeline uses **DESeq2** to identify differentially expressed genes across conditions and outputs various visualizations.
- **Co-expression Network Analysis (WGCNA)**: **WGCNA** is used to identify modules of co-expressed genes and generate a network structure based on expression similarities.
- **Interactive Visualizations**: Generates easy-to-interpret visual outputs, including PCA, heatmaps, volcano plots, and MA plots.
- **Gene ID Mapping**: Automatically maps **NCBI Entrez Gene IDs** to **HGNC gene symbols** using the **org.Hs.eg.db** database for easy interpretation.
---

#### Example Outputs:

##### 1. **Volcano Plot**: Visualizing differentially expressed genes

<p align="center">
  <img src="volcano_plot_tumor_vs_normal.png" alt="Volcano Plot" width="600"/>
</p>

The volcano plot shows significant genes based on log2 fold change and adjusted p-values (padj). Genes with higher fold changes and statistical significance are highlighted in red.

##### 2. **Principal Component Analysis (PCA) Plot**: Visualizing overall variance between samples

<p align="center">
  <img src="pca_plot.png" alt="PCA Plot" width="600"/>
</p>

PCA plot of the variance-stabilized data helps to assess sample separation based on the top principal components, offering insight into the major sources of variance in the data.

##### 3. **MA Plot**: Log fold change vs. mean expression

<p align="center">
  <img src="MA_plot.png" alt="MA Plot" width="600"/>
</p>

The MA plot visualizes log2 fold changes against the mean expression for all genes, where significant genes are highlighted.
##### 4. **Heatmap**: Visualization of Top Differentially Expressed Genes

<p align="center">
  <img src="heatmap.png" alt="Heatmap" width="600"/>
</p>

The heatmap provides a visualization of the expression patterns of the top 50 differentially expressed genes across samples. Each row represents a gene, and each column represents a sample. The color intensity indicates the relative expression level after centering (subtracting the row mean). Samples are clustered by similarity to reveal patterns associated with different conditions, such as “normal” and “tumor,” allowing for quick identification of co-expression patterns and differences between groups.

##### 5. **Co-expression Network (WGCNA) Dendrogram**: Module detection via hierarchical clustering

<p align="center">
  <img src="dendrogram_with_modules.png" alt="WGCNA Dendrogram" width="600"/>
</p>

This dendrogram displays co-expressed gene modules identified through hierarchical clustering in WGCNA. Modules are color-coded based on gene similarity, providing a visual grouping of genes that are likely co-regulated. Each branch represents a module with genes sharing similar expression patterns, potentially associated with specific biological functions.

##### 6. **Module-Trait Relationship Heatmap**: Correlating Modules with External Traits

<p align="center">
  <img src="module_trait_relationships_conditions_no_labels_adjusted.png" alt="Module-Trait Relationship Heatmap
" width="600"/>
</p>

The module-trait relationship heatmap provides a summary of the correlations between gene modules and external sample traits (e.g., tumor vs. normal). Each row represents a module (colored by module name, such as "blue" or "turquoise"), and each column represents a sample trait. The color of each cell reflects the strength and direction of the correlation: red for positive correlations, blue for negative correlations, and white for no correlation.

The heatmap legend indicates the correlation scale from -1 to 1. Modules with strong positive or negative correlations to a trait may contain genes involved in processes related to that condition. This heatmap allows for quick identification of potentially biologically relevant modules associated with specific sample traits, making it a useful tool for further investigation into gene-trait relationships.

---

#### Workflow:

1. **Data Preprocessing**:
   - Download sample metadata from GEO.
   - Load the RNA-seq count data.
   - Map gene IDs (Entrez to HGNC symbols).

2. **Differential Expression**:
   - Run **DESeq2** to perform differential expression analysis between conditions.
   - Output results as CSV and generate a volcano plot to visualize significant genes.

3. **Normalization**:
   - Normalize the RNA-seq data using **variance stabilizing transformation (VST)**.

4. **Co-expression Analysis**:
   - Perform **WGCNA** to detect gene modules that are co-expressed.
   - Output dendrogram plots and visualize co-expression networks.
---

#### Tools and Technologies:
- **R**: The pipeline is built using R and Bioconductor packages.
- **Packages Used**:
  - **DESeq2**: For differential expression analysis and normalization.
  - **WGCNA**: For co-expression network analysis.
  - **GEOquery**: To download and process GEO datasets.
  - **ggplot2**, **pheatmap**: For generating heatmaps and visualizations.
  - **EnhancedVolcano**: For creating volcano plots.
  - **clusterProfiler**: For KEGG pathway enrichment analysis.
  - **org.Hs.eg.db**: For gene ID mapping.

---

#### Results & Deliverables:
- **differential_expression_results.csv**: Contains log fold changes, p-values, and adjusted p-values for each gene.
- **Volcano Plot, PCA Plot, Heatmap, MA Plot**: High-quality figures summarizing the analysis.
- **Co-expression Network Results**: Identifies gene modules and network structure.
---

#### Example Usage:

```r
# Example Usage for GEO dataset GSE37764
geo_id <- "GSE37764"
counts_file <- "path_to/GSE37764_raw_counts_GRCh38.p13_NCBI.tsv"
```

---

