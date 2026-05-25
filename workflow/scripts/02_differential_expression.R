# 02_differential_expression.R
# DESeq2 differential expression for treatment vs reference, with lfcShrink,
# a VST object cached for WGCNA, and the standard volcano / MA / heatmap.

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
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
# Don't trust on-disk row order to match counts column order — reorder + verify.
samples <- samples[match(colnames(counts), samples$sample_id), , drop = FALSE]
stopifnot(!anyNA(samples$sample_id),
          identical(samples$sample_id, colnames(counts)))

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
# apeglm requires `coef` to exist verbatim in resultsNames(dds). DESeq2 mangles
# level names (spaces, punctuation) when it builds coefficient names, so verify
# rather than trust string interpolation.
coef_name <- paste0("condition_", cfg$dataset$treatment_level,
                    "_vs_", cfg$dataset$reference_level)
if (!coef_name %in% resultsNames(dds)) {
  stop("lfcShrink coef '", coef_name, "' not found in resultsNames(dds): ",
       paste(resultsNames(dds), collapse = ", "))
}
res <- lfcShrink(dds, coef = coef_name, res = res,
                 type = cfg$deseq2$shrink_method)
res <- res[order(res$padj), ]

# Attach gene symbols from genes.tsv (written by 01) instead of re-querying
# org.Hs.eg.db — recount3's rowData already carries authoritative GENCODE
# gene_name values keyed on the same versioned ENSEMBL IDs.
genes_path <- file.path(dirname(cfg$paths$raw_counts), "genes.tsv")
if (file.exists(genes_path)) {
  genes <- read.table(genes_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
  res$gene_symbol <- genes$gene_name[match(rownames(res), genes$gene_id)]
} else {
  res$gene_symbol <- NA_character_
}

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
  breaks = seq(-3, 3, length.out = 101),   # clip extremes so a few cells don't wash out the palette
  main = paste0("Top ", cfg$deseq2$top_n_heatmap,
                " DE genes — ", cfg$project$cancer_short),
  filename = file.path(cfg$paths$figures_dir, "heatmap_top_de.png"),
  width = 10, height = 10
)
