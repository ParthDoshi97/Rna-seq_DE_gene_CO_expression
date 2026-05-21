# 02_differential_expression.R
# DESeq2 differential expression for treatment vs reference, with lfcShrink,
# a VST object cached for WGCNA, and the standard volcano / MA / heatmap.

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
cfg <- yaml::read_yaml(config_path)

set.seed(cfg$runtime$seed)
dir.create(cfg$paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$rds_dir,     recursive = TRUE, showWarnings = FALSE)

# ---- Load counts + samples --------------------------------------------------

counts_df <- read.table(cfg$paths$raw_counts, header = TRUE, sep = "\t",
                        check.names = FALSE)
rownames(counts_df) <- counts_df$gene_id
counts <- as.matrix(counts_df[, -1, drop = FALSE])

samples <- read.table(cfg$paths$sample_metadata, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
stopifnot(all(samples$sample_id == colnames(counts)))

samples$condition <- factor(
  samples$condition,
  levels = c(cfg$dataset$reference_level, cfg$dataset$treatment_level)
)

# ---- DESeq2 -----------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = samples,
                              design    = ~ condition)
dds <- DESeq(dds, parallel = FALSE)

res <- results(
  dds,
  contrast = c("condition", cfg$dataset$treatment_level, cfg$dataset$reference_level),
  alpha = cfg$deseq2$padj_threshold,
  independentFiltering = cfg$deseq2$independent_filtering
)
res <- lfcShrink(dds,
                 coef = paste0("condition_", cfg$dataset$treatment_level,
                               "_vs_", cfg$dataset$reference_level),
                 res = res, type = cfg$deseq2$shrink_method)
res <- res[order(res$padj), ]

# Map gene IDs → symbols. recount3 IDs are ENSEMBL (versioned).
ensembl_ids <- sub("\\..*$", "", rownames(res))
res$gene_symbol <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                          column = "SYMBOL", keytype = "ENSEMBL",
                          multiVals = "first")

# Save dds + vsd + res for downstream
vsd <- vst(dds, blind = FALSE)
saveRDS(dds, cfg$paths$dds_rds)
saveRDS(vsd, cfg$paths$vsd_rds)
saveRDS(res, cfg$paths$de_results_rds)

# ---- Deliverable table ------------------------------------------------------

de_out <- data.frame(
  gene_id     = rownames(res),
  gene_symbol = res$gene_symbol,
  baseMean    = res$baseMean,
  log2FoldChange = res$log2FoldChange,
  lfcSE       = res$lfcSE,
  pvalue      = res$pvalue,
  padj        = res$padj,
  stringsAsFactors = FALSE
)
write.table(de_out, file = cfg$paths$de_results_csv,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

n_sig_up   <- sum(de_out$padj < cfg$deseq2$padj_threshold &
                  de_out$log2FoldChange >  cfg$deseq2$log2fc_threshold, na.rm = TRUE)
n_sig_down <- sum(de_out$padj < cfg$deseq2$padj_threshold &
                  de_out$log2FoldChange < -cfg$deseq2$log2fc_threshold, na.rm = TRUE)
message("DE genes (padj < ", cfg$deseq2$padj_threshold,
        ", |log2FC| > ", cfg$deseq2$log2fc_threshold, "): ",
        n_sig_up, " up, ", n_sig_down, " down")

# ---- Figures ----------------------------------------------------------------

# Volcano
volcano <- EnhancedVolcano(
  res,
  lab = res$gene_symbol,
  x = "log2FoldChange", y = "padj",
  pCutoff  = cfg$deseq2$padj_threshold,
  FCcutoff = cfg$deseq2$log2fc_threshold,
  pointSize = 2.0, labSize = 3.0,
  title = paste0("DE: ", cfg$dataset$treatment_level, " vs ",
                 cfg$dataset$reference_level, " — ", cfg$project$cancer_short),
  subtitle = "DESeq2 (apeglm shrink)"
)
ggsave(file.path(cfg$paths$figures_dir, "volcano.png"),
       volcano, width = 8, height = 6, dpi = 150)

# MA
png(file.path(cfg$paths$figures_dir, "ma_plot.png"),
    width = 8, height = 6, units = "in", res = 150)
DESeq2::plotMA(res, ylim = c(-5, 5),
               main = paste0("MA — ", cfg$project$cancer_short))
dev.off()

# Heatmap of top-N DE genes
top_idx <- head(order(res$padj), cfg$deseq2$top_n_heatmap)
heat_mat <- assay(vsd)[top_idx, , drop = FALSE]
heat_mat <- heat_mat - rowMeans(heat_mat)
rownames(heat_mat) <- ifelse(is.na(res$gene_symbol[top_idx]),
                             rownames(res)[top_idx],
                             res$gene_symbol[top_idx])
pheatmap(
  heat_mat,
  annotation_col = data.frame(condition = samples$condition,
                              row.names = samples$sample_id),
  show_rownames = TRUE, show_colnames = FALSE,
  main = paste0("Top ", cfg$deseq2$top_n_heatmap,
                " DE genes — ", cfg$project$cancer_short),
  filename = file.path(cfg$paths$figures_dir, "heatmap_top_de.png"),
  width = 10, height = 10
)
