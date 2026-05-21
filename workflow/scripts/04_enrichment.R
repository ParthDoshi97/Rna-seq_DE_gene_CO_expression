#!/usr/bin/env Rscript
# 04_enrichment.R
# clusterProfiler GO (BP) + KEGG enrichment for up- and down-regulated
# DE genes. Per-module enrichment is added in branch 5.

suppressPackageStartupMessages({
  library(yaml)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(enrichplot)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
cfg <- yaml::read_yaml(config_path)

set.seed(cfg$runtime$seed)
dir.create(cfg$paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$tables_dir,  recursive = TRUE, showWarnings = FALSE)

# ---- Inputs -----------------------------------------------------------------

res <- readRDS(cfg$paths$de_results_rds)
de <- as.data.frame(res)
de$gene_id <- rownames(res)
de$ensembl <- sub("\\..*$", "", de$gene_id)
de$entrez  <- mapIds(org.Hs.eg.db, keys = de$ensembl,
                     column = "ENTREZID", keytype = "ENSEMBL",
                     multiVals = "first")

universe <- unique(na.omit(de$entrez))

sig <- de$padj < cfg$deseq2$padj_threshold & !is.na(de$padj)
up_genes   <- unique(na.omit(de$entrez[sig & de$log2FoldChange >  cfg$deseq2$log2fc_threshold]))
down_genes <- unique(na.omit(de$entrez[sig & de$log2FoldChange < -cfg$deseq2$log2fc_threshold]))
message("Enrichment inputs: ", length(up_genes), " up, ",
        length(down_genes), " down (universe: ", length(universe), ")")

# ---- Helper -----------------------------------------------------------------

run_go <- function(genes, label) {
  if (length(genes) < cfg$enrichment$min_gs_size) return(NULL)
  enrichGO(gene = genes, universe = universe,
           OrgDb = cfg$enrichment$org_db,
           keyType = cfg$enrichment$keytype,
           ont = cfg$enrichment$go_ontology,
           pvalueCutoff  = cfg$enrichment$pvalue_cutoff,
           qvalueCutoff  = cfg$enrichment$qvalue_cutoff,
           minGSSize = cfg$enrichment$min_gs_size,
           maxGSSize = cfg$enrichment$max_gs_size,
           readable = TRUE)
}

run_kegg <- function(genes, label) {
  genes <- unique(as.character(genes))
  bg <- unique(as.character(universe))
  
  if (length(genes) < cfg$enrichment$min_gs_size) {
    message("  skipping KEGG for ", label, ": too few genes")
    return(NULL)
  }
  
  tryCatch({
    clusterProfiler::enrichKEGG(
      gene = genes,
      universe = bg,
      organism = cfg$enrichment$kegg_organism,
      keyType = "ncbi-geneid",
      pvalueCutoff = cfg$enrichment$pvalue_cutoff,
      qvalueCutoff = cfg$enrichment$qvalue_cutoff,
      minGSSize = cfg$enrichment$min_gs_size,
      maxGSSize = cfg$enrichment$max_gs_size
    )
  }, error = function(e) {
    message("  KEGG failed for ", label, ": ", conditionMessage(e))
    return(NULL)
  })
}

save_dotplot <- function(obj, file, title) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    message("  no enriched terms for ", title)
    return(invisible())
  }
  p <- dotplot(obj, showCategory = cfg$enrichment$top_n_dotplot,
               title = title) +
    theme(plot.title = element_text(size = 11))
  ggsave(file, p, width = 9, height = 7, dpi = 150)
}

# ---- Run --------------------------------------------------------------------

go_up   <- run_go(up_genes,   "GO BP — up")
go_down <- run_go(down_genes, "GO BP — down")
kegg_up   <- run_kegg(up_genes,   "KEGG — up")
kegg_down <- run_kegg(down_genes, "KEGG — down")

save_dotplot(go_up,   file.path(cfg$paths$figures_dir, "go_bp_up.png"),
             "GO BP enrichment — up-regulated")
save_dotplot(go_down, file.path(cfg$paths$figures_dir, "go_bp_down.png"),
             "GO BP enrichment — down-regulated")
save_dotplot(kegg_up,   file.path(cfg$paths$figures_dir, "kegg_up.png"),
             "KEGG enrichment — up-regulated")
save_dotplot(kegg_down, file.path(cfg$paths$figures_dir, "kegg_down.png"),
             "KEGG enrichment — down-regulated")

# ---- Tables -----------------------------------------------------------------

bind_rows_safe <- function(...) {
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, list(...))
  if (length(parts) == 0) return(data.frame())
  do.call(rbind, parts)
}

go_tbl <- bind_rows_safe(
  if (!is.null(go_up))   cbind(direction = "up",   as.data.frame(go_up)),
  if (!is.null(go_down)) cbind(direction = "down", as.data.frame(go_down))
)
kegg_tbl <- bind_rows_safe(
  if (!is.null(kegg_up))   cbind(direction = "up",   as.data.frame(kegg_up)),
  if (!is.null(kegg_down)) cbind(direction = "down", as.data.frame(kegg_down))
)

write.table(go_tbl,   file = cfg$paths$enrichment_go_csv,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(kegg_tbl, file = cfg$paths$enrichment_kegg_csv,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
