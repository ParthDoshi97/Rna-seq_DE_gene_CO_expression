#!/usr/bin/env Rscript
# 04_enrichment.R
# clusterProfiler functional enrichment for tumor-vs-normal DE results.
#   - ORA  (over-representation) on up/down gene lists  -- GO + KEGG
#   - GSEA on the LFC-ranked gene list                  -- GO + KEGG
# Per-module enrichment is added in branch 5.

suppressPackageStartupMessages({
  library(yaml)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(enrichplot)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
cfg <- yaml::read_yaml(config_path)

# enrichGO's `keyType` must match the type of de$entrez (Entrez IDs). A stale
# config that still said "ENSEMBL" would silently fail to map every gene and
# return an empty result table for no obvious reason.
stopifnot(cfg$enrichment$keytype == "ENTREZID")

set.seed(cfg$runtime$seed)
dir.create(cfg$paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$tables_dir,  recursive = TRUE, showWarnings = FALSE)

# ---- Inputs -----------------------------------------------------------------

res <- readRDS(cfg$paths$de_results_rds)
de  <- as.data.frame(res)
de$gene_id <- rownames(res)
de$ensembl <- sub("\\..*$", "", de$gene_id)

# multiVals = "first" picks one Entrez ID per Ensembl when the mapping is 1:N.
# Fine for enrichment-level analysis, but it's a lossy choice -- runner-up
# isoforms are dropped.
de$entrez <- mapIds(org.Hs.eg.db, keys = de$ensembl,
                    column = "ENTREZID", keytype = "ENSEMBL",
                    multiVals = "first")

n_de     <- nrow(de)
n_mapped <- sum(!is.na(de$entrez))
message("ENSEMBL->ENTREZ: ", n_mapped, " / ", n_de,
        " mapped (", round(100 * n_mapped / max(n_de, 1), 1), "%)")

universe <- unique(na.omit(de$entrez))

sig <- de$padj < cfg$deseq2$padj_threshold & !is.na(de$padj)
up_genes   <- unique(na.omit(de$entrez[sig & de$log2FoldChange >  cfg$deseq2$log2fc_threshold]))
down_genes <- unique(na.omit(de$entrez[sig & de$log2FoldChange < -cfg$deseq2$log2fc_threshold]))
message("Enrichment inputs: ", length(up_genes), " up, ",
        length(down_genes), " down (universe: ", length(universe), ")")

# ---- GO ontologies + thresholds --------------------------------------------
# min_gs_size and min_input_genes are conceptually different:
#   min_gs_size      lower bound on a *gene set's* member count (clusterProfiler arg)
#   min_input_genes  lower bound on the *input gene list* -- below this we skip the run

min_input     <- cfg$enrichment$min_input_genes %||% 10
go_ontologies <- cfg$enrichment$go_ontology
if (length(go_ontologies) == 1 && go_ontologies == "ALL") {
  go_ontologies <- c("BP", "MF", "CC")
}

# ---- ORA helpers ------------------------------------------------------------

run_go <- function(genes, label, ont) {
  if (length(genes) < min_input) {
    message("  skipping GO ", ont, " (", label, "): ",
            length(genes), " < min_input_genes=", min_input)
    return(NULL)
  }
  message("  GO ", ont, " (", label, "): ", length(genes), " gene(s)")
  enrichGO(gene = genes, universe = universe,
           OrgDb = cfg$enrichment$org_db,
           keyType = cfg$enrichment$keytype,
           ont = ont,
           pvalueCutoff = cfg$enrichment$pvalue_cutoff,
           qvalueCutoff = cfg$enrichment$qvalue_cutoff,
           minGSSize = cfg$enrichment$min_gs_size,
           maxGSSize = cfg$enrichment$max_gs_size,
           readable = TRUE)
}

run_kegg <- function(genes, label) {
  genes <- unique(as.character(genes))
  bg    <- unique(as.character(universe))
  if (length(genes) < min_input) {
    message("  skipping KEGG (", label, "): ",
            length(genes), " < min_input_genes=", min_input)
    return(NULL)
  }
  message("  KEGG (", label, "): ", length(genes), " gene(s)")
  tryCatch({
    clusterProfiler::enrichKEGG(
      gene = genes, universe = bg,
      organism = cfg$enrichment$kegg_organism,
      keyType = "ncbi-geneid",
      pvalueCutoff = cfg$enrichment$pvalue_cutoff,
      qvalueCutoff = cfg$enrichment$qvalue_cutoff,
      minGSSize = cfg$enrichment$min_gs_size,
      maxGSSize = cfg$enrichment$max_gs_size
    )
  }, error = function(e) {
    message("  KEGG failed (", label, "): ", conditionMessage(e))
    NULL
  })
}

# ---- GSEA helpers -----------------------------------------------------------
# ORA collapses everything to "above / below a hard cutoff". With a cohort
# this size most genes are "significant" and ORA recovers only generic terms
# ("cell cycle", "extracellular matrix"). GSEA on the LFC-ranked list uses
# the whole gradient and usually resolves cleaner biology.

build_ranked <- function() {
  ok <- !is.na(de$entrez) & !is.na(de$log2FoldChange)
  v  <- de$log2FoldChange[ok]
  names(v) <- de$entrez[ok]
  # Multiple Ensembl IDs collapsing to the same Entrez: keep the strongest |LFC|.
  v <- tapply(v, names(v), function(x) x[which.max(abs(x))])
  sort(v, decreasing = TRUE)
}

run_gsea_go <- function(ranked, ont) {
  message("  GSEA GO ", ont, ": ", length(ranked), " ranked gene(s)")
  tryCatch({
    gseGO(geneList = ranked,
          OrgDb    = cfg$enrichment$org_db,
          keyType  = cfg$enrichment$keytype,
          ont      = ont,
          minGSSize    = cfg$enrichment$min_gs_size,
          maxGSSize    = cfg$enrichment$max_gs_size,
          pvalueCutoff = cfg$enrichment$pvalue_cutoff,
          verbose      = FALSE)
  }, error = function(e) {
    message("  GSEA GO ", ont, " failed: ", conditionMessage(e))
    NULL
  })
}

run_gsea_kegg <- function(ranked) {
  message("  GSEA KEGG: ", length(ranked), " ranked gene(s)")
  tryCatch({
    gseKEGG(geneList = ranked,
            organism = cfg$enrichment$kegg_organism,
            keyType  = "ncbi-geneid",
            minGSSize    = cfg$enrichment$min_gs_size,
            maxGSSize    = cfg$enrichment$max_gs_size,
            pvalueCutoff = cfg$enrichment$pvalue_cutoff,
            verbose      = FALSE)
  }, error = function(e) {
    message("  GSEA KEGG failed: ", conditionMessage(e))
    NULL
  })
}

# ---- Dotplot helper ---------------------------------------------------------

save_dotplot <- function(obj, file, title, gsea = FALSE) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    message("  no enriched terms for ", title)
    return(invisible())
  }
  p <- if (gsea) {
    dotplot(obj, showCategory = cfg$enrichment$top_n_dotplot,
            title = title, split = ".sign") + facet_grid(. ~ .sign)
  } else {
    dotplot(obj, showCategory = cfg$enrichment$top_n_dotplot, title = title)
  }
  p <- p + theme(plot.title = element_text(size = 11))
  ggsave(file, p, width = if (gsea) 11 else 9, height = 7, dpi = 150)
}

# ---- Run ORA ---------------------------------------------------------------

message("Running ORA...")
go_up   <- lapply(go_ontologies, function(o) run_go(up_genes,   "up",   o))
go_down <- lapply(go_ontologies, function(o) run_go(down_genes, "down", o))
names(go_up) <- names(go_down) <- go_ontologies
kegg_up   <- run_kegg(up_genes,   "up")
kegg_down <- run_kegg(down_genes, "down")

# ---- Run GSEA --------------------------------------------------------------

message("Running GSEA...")
ranked    <- build_ranked()
gsea_go   <- lapply(go_ontologies, function(o) run_gsea_go(ranked, o))
names(gsea_go) <- go_ontologies
gsea_kegg <- run_gsea_kegg(ranked)

# ---- Reproducibility breadcrumb for KEGG -----------------------------------
# clusterProfiler hits KEGG REST at runtime, so results drift as KEGG updates.
# A sidecar with package versions + timestamp keeps the run traceable in the
# audit trail without needing internet on re-read.

kegg_meta_path <- file.path(cfg$paths$tables_dir, "kegg_run_info.txt")
writeLines(c(
  paste0("kegg_organism: ", cfg$enrichment$kegg_organism),
  paste0("clusterProfiler_version: ",
         as.character(packageVersion("clusterProfiler"))),
  paste0("org_db_version: ",
         as.character(packageVersion(cfg$enrichment$org_db))),
  paste0("queried_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R_version: ", R.version.string)
), kegg_meta_path)
message("Wrote ", kegg_meta_path)

# ---- Figures ---------------------------------------------------------------

for (o in go_ontologies) {
  save_dotplot(go_up[[o]],
               file.path(cfg$paths$figures_dir, paste0("go_", tolower(o), "_up.png")),
               paste0("GO ", o, " enrichment - up-regulated"))
  save_dotplot(go_down[[o]],
               file.path(cfg$paths$figures_dir, paste0("go_", tolower(o), "_down.png")),
               paste0("GO ", o, " enrichment - down-regulated"))
  save_dotplot(gsea_go[[o]],
               file.path(cfg$paths$figures_dir, paste0("gsea_go_", tolower(o), ".png")),
               paste0("GSEA GO ", o), gsea = TRUE)
}
save_dotplot(kegg_up,
             file.path(cfg$paths$figures_dir, "kegg_up.png"),
             "KEGG enrichment - up-regulated")
save_dotplot(kegg_down,
             file.path(cfg$paths$figures_dir, "kegg_down.png"),
             "KEGG enrichment - down-regulated")
save_dotplot(gsea_kegg,
             file.path(cfg$paths$figures_dir, "gsea_kegg.png"),
             "GSEA KEGG", gsea = TRUE)

# ---- Tables ----------------------------------------------------------------

# Different ontologies/directions can return different optional columns; align
# by union so rbind doesn't drop any of them.
bind_rows_safe <- function(parts) {
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, parts)
  if (length(parts) == 0) return(data.frame())
  all_cols <- unique(unlist(lapply(parts, colnames)))
  parts <- lapply(parts, function(x) {
    miss <- setdiff(all_cols, colnames(x))
    for (m in miss) x[[m]] <- NA
    x[, all_cols, drop = FALSE]
  })
  do.call(rbind, parts)
}

go_tbl_parts <- list()
for (o in go_ontologies) {
  if (!is.null(go_up[[o]]))
    go_tbl_parts[[paste0(o, "_up")]]   <- cbind(direction = "up",   ontology = o,
                                                as.data.frame(go_up[[o]]))
  if (!is.null(go_down[[o]]))
    go_tbl_parts[[paste0(o, "_down")]] <- cbind(direction = "down", ontology = o,
                                                as.data.frame(go_down[[o]]))
}
go_tbl <- bind_rows_safe(go_tbl_parts)

kegg_tbl <- bind_rows_safe(list(
  if (!is.null(kegg_up))   cbind(direction = "up",   as.data.frame(kegg_up)),
  if (!is.null(kegg_down)) cbind(direction = "down", as.data.frame(kegg_down))
))

gsea_go_parts <- list()
for (o in go_ontologies) {
  if (!is.null(gsea_go[[o]]))
    gsea_go_parts[[o]] <- cbind(ontology = o, as.data.frame(gsea_go[[o]]))
}
gsea_go_tbl   <- bind_rows_safe(gsea_go_parts)
gsea_kegg_tbl <- if (!is.null(gsea_kegg)) as.data.frame(gsea_kegg) else data.frame()

write.table(go_tbl,   file = cfg$paths$enrichment_go_csv,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(kegg_tbl, file = cfg$paths$enrichment_kegg_csv,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

gsea_go_path   <- cfg$paths$gsea_go_csv   %||% file.path(cfg$paths$tables_dir, "gsea_go.tsv")
gsea_kegg_path <- cfg$paths$gsea_kegg_csv %||% file.path(cfg$paths$tables_dir, "gsea_kegg.tsv")
write.table(gsea_go_tbl,   file = gsea_go_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(gsea_kegg_tbl, file = gsea_kegg_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

message("Done.")
