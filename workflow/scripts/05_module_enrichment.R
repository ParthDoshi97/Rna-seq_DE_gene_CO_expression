# 05_module_enrichment.R
# Functionally annotate each WGCNA module via GO/KEGG over-representation
# analysis. Tests modules against the background universe of variance-selected
# genes (the network's own input: ~5000 genes), NOT the whole genome.

suppressPackageStartupMessages({
  library(yaml)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
cfg <- yaml::read_yaml(config_path)

# All paths in config.yaml are relative to the repo root.
# Resolve repo root from config file location and setwd so relative paths work
# regardless of whether the script is run from workflow/scripts/ or repo root.
config_abs <- normalizePath(config_path, mustWork = TRUE)
repo_root  <- dirname(dirname(config_abs))   # <repo>/config/config.yaml -> <repo>
setwd(repo_root)
message("Working directory set to repo root: ", repo_root)

set.seed(cfg$runtime$seed)
dir.create(cfg$paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$rds_dir,     recursive = TRUE, showWarnings = FALSE)

# ---- Load inputs ----

wgcna <- readRDS(cfg$paths$wgcna_rds)
moduleColors <- wgcna$moduleColors
datExpr <- wgcna$datExpr

module_assignments_path <- cfg$paths$module_assignments_csv %||%
  file.path(cfg$paths$tables_dir, "module_assignments.tsv")
module_assignments <- read.table(module_assignments_path,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

genes_path <- file.path(dirname(cfg$paths$raw_counts), "genes.tsv")
if (!file.exists(genes_path)) {
  stop("genes.tsv not found at ", genes_path)
}
genes <- read.table(genes_path, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)

# Config
me_cfg <- cfg$module_enrichment %||% list()
keytype <- me_cfg$keytype %||% "ENTREZID"
go_ontologies <- me_cfg$go_ontologies %||% "BP"
pvalue_cutoff <- me_cfg$pvalue_cutoff %||% 0.05
qvalue_cutoff <- me_cfg$qvalue_cutoff %||% 0.10
min_module_size <- me_cfg$min_module_size %||% 30
min_gs_size <- me_cfg$min_gs_size %||% 10
max_gs_size <- me_cfg$max_gs_size %||% 500
run_kegg <- me_cfg$run_kegg %||% TRUE
top_terms_plot <- me_cfg$top_terms_plot %||% 15
label_wrap_width <- me_cfg$label_wrap_width %||% 34
plot_width <- me_cfg$plot_width %||% 12
plot_min_height <- me_cfg$plot_min_height %||% 6
plot_row_height <- me_cfg$plot_row_height %||% 0.42
plot_extra_line_height <- me_cfg$plot_extra_line_height %||% 0.18

stopifnot("keytype must be ENTREZID" = keytype == "ENTREZID")

wrap_term_labels <- function(labels, width) {
  vapply(labels, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  }, character(1))
}

plot_height_for_labels <- function(wrapped_labels) {
  n_lines <- vapply(strsplit(wrapped_labels, "\n", fixed = TRUE),
                    length, integer(1))
  extra_lines <- sum(pmax(n_lines - 1, 0))
  max(plot_min_height,
      2.6 + plot_row_height * length(wrapped_labels) +
        plot_extra_line_height * extra_lines)
}

# ---- Build background universe ----

universe_genes <- unique(module_assignments$gene_id)
message("Universe: ", length(universe_genes), " genes from variance-selected input")

# Strip version numbers from ENSEMBL IDs (e.g., ENSG00000187634.11 -> ENSG00000187634)
universe_genes_clean <- sub("\\..*$", "", universe_genes)

universe_entrez <- mapIds(org.Hs.eg.db,
                          keys = universe_genes_clean,
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")
universe_entrez <- na.omit(universe_entrez)
universe_entrez <- as.character(universe_entrez)

message("Universe mapped to ENTREZ: ", length(universe_entrez), " / ",
        length(universe_genes), " (", round(100 * length(universe_entrez) / length(universe_genes), 1), "%)")
stopifnot("Failed to map universe to ENTREZ" = length(universe_entrez) > 0)

# ---- Build per-module gene lists ----

modules <- setdiff(unique(module_assignments$module), "grey")
module_genes <- split(module_assignments$gene_id,
                      module_assignments$module)
module_genes <- module_genes[modules]

message("Modules: ", length(modules), " (excluding grey)")

module_entrez <- list()
for (m in modules) {
  genes_in_m <- module_genes[[m]]
  # Strip version numbers from ENSEMBL IDs
  genes_in_m_clean <- sub("\\..*$", "", genes_in_m)
  entrez <- mapIds(org.Hs.eg.db,
                   keys = genes_in_m_clean,
                   column = "ENTREZID",
                   keytype = "ENSEMBL",
                   multiVals = "first")
  entrez <- na.omit(entrez)
  entrez <- as.character(entrez)

  if (length(entrez) >= min_module_size) {
    module_entrez[[m]] <- entrez
  } else {
    message("  Skipped ", m, ": ", length(entrez), " genes < min_module_size (",
            min_module_size, ")")
  }
}

message("Modules to test: ", length(module_entrez))

# ---- Run enrichment per module ----

go_results <- list()
kegg_results <- list()

for (m in names(module_entrez)) {
  genes_m <- module_entrez[[m]]

  for (ont in go_ontologies) {
    tryCatch({
      ego <- enrichGO(gene = genes_m,
                      universe = universe_entrez,
                      OrgDb = org.Hs.eg.db,
                      keyType = keytype,
                      ont = ont,
                      pvalueCutoff = pvalue_cutoff,
                      qvalueCutoff = qvalue_cutoff,
                      minGSSize = min_gs_size,
                      maxGSSize = max_gs_size,
                      readable = TRUE)
      if (!is.null(ego) && nrow(ego) > 0) {
        ego_df <- as.data.frame(ego)
        ego_df$module <- m
        ego_df$ontology <- ont
        go_results[[paste0(m, "_", ont)]] <- ego_df
      }
    }, error = function(e) {
      warning("GO enrichment failed for module ", m, " / ", ont, ": ",
              e$message)
    })
  }

  if (run_kegg) {
    tryCatch({
      ekegg <- enrichKEGG(gene = genes_m,
                          universe = universe_entrez,
                          organism = "hsa",
                          pvalueCutoff = pvalue_cutoff,
                          minGSSize = min_gs_size,
                          maxGSSize = max_gs_size)
      if (!is.null(ekegg) && nrow(ekegg) > 0) {
        ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = keytype)
        ekegg_df <- as.data.frame(ekegg)
        ekegg_df$module <- m
        kegg_results[[m]] <- ekegg_df
      }
    }, error = function(e) {
      warning("KEGG enrichment failed for module ", m, ": ", e$message)
    })
  }
}

# ---- Assemble results ----

if (length(go_results) > 0) {
  go_tbl <- do.call(rbind, go_results)
  rownames(go_tbl) <- NULL
} else {
  go_tbl <- data.frame(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer(),
    module = character(), ontology = character()
  )
}

if (length(kegg_results) > 0) {
  kegg_tbl <- do.call(rbind, kegg_results)
  rownames(kegg_tbl) <- NULL
} else {
  kegg_tbl <- data.frame(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer(),
    module = character()
  )
}

# ---- Module summary table ----

summary_rows <- list()
for (m in names(module_entrez)) {
  m_go <- go_tbl[go_tbl$module == m & !is.na(go_tbl$ID), ]
  m_go <- m_go[order(m_go$p.adjust), ]

  size <- length(module_entrez[[m]])
  top_terms <- if (nrow(m_go) > 0) {
    head(m_go$Description, 3)
  } else {
    rep(NA_character_, 3)
  }
  best_padj <- if (nrow(m_go) > 0) m_go$p.adjust[1] else NA_real_

  summary_rows[[m]] <- data.frame(
    module = m,
    size = size,
    top_term_1 = top_terms[1],
    top_term_2 = if (length(top_terms) > 1) top_terms[2] else NA_character_,
    top_term_3 = if (length(top_terms) > 2) top_terms[3] else NA_character_,
    best_padj = best_padj,
    stringsAsFactors = FALSE
  )
}
summary_tbl <- do.call(rbind, summary_rows)
rownames(summary_tbl) <- NULL

# ---- Write tables ----

go_path <- cfg$paths$module_enrichment_go_csv %||%
  file.path(cfg$paths$tables_dir, "module_enrichment_go.tsv")
kegg_path <- cfg$paths$module_enrichment_kegg_csv %||%
  file.path(cfg$paths$tables_dir, "module_enrichment_kegg.tsv")
summary_path <- cfg$paths$module_annotation_summary_csv %||%
  file.path(cfg$paths$tables_dir, "module_annotation_summary.tsv")

write.table(go_tbl, file = go_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(kegg_tbl, file = kegg_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(summary_tbl, file = summary_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

message("Wrote GO results: ", go_path)
message("Wrote KEGG results: ", kegg_path)
message("Wrote summary: ", summary_path)

# ---- Figures ----

fig_dir <- file.path(cfg$paths$figures_dir, "05_module_enrichment")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

for (m in names(module_entrez)) {
  m_go <- go_tbl[go_tbl$module == m & !is.na(go_tbl$ID), ]

  if (nrow(m_go) == 0) {
    message("  Skipped figure for ", m, ": no significant terms")
    next
  }

  m_go <- m_go[order(m_go$p.adjust), ]
  m_go <- head(m_go, top_terms_plot)

  tryCatch({
    m_go$Description <- factor(m_go$Description,
                               levels = rev(m_go$Description))
    wrapped_labels <- wrap_term_labels(as.character(m_go$Description),
                                       label_wrap_width)
    plot_height <- plot_height_for_labels(wrapped_labels)
    png(file.path(fig_dir, paste0("enrichment_", m, ".png")),
        width = plot_width, height = plot_height,
        units = "in", res = 150)
    par(mar = c(5, 18, 4, 2))
    barplot(rev(-log10(m_go$p.adjust)), names.arg = rev(wrapped_labels),
            horiz = TRUE, las = 1, col = m, border = NA,
            xlab = "-log10(p.adjust)",
            main = paste0("Module enrichment - ", m),
            cex.names = 0.78, cex.axis = 0.9,
            cex.lab = 0.95, cex.main = 1.05)
    dev.off()
    message("Wrote figure: ", fig_dir, "/enrichment_", m, ".png")
  }, error = function(e) {
    warning("Figure generation failed for ", m, ": ", e$message)
  })
}

# ---- Write run info ----

run_info <- paste0(
  "Module Enrichment Run Info\n",
  "==========================\n",
  "Timestamp: ", Sys.time(), "\n",
  "Script: 05_module_enrichment.R\n",
  "clusterProfiler version: ", packageVersion("clusterProfiler"), "\n",
  "org.Hs.eg.db version: ", packageVersion("org.Hs.eg.db"), "\n",
  "\n",
  "Universe: ", length(universe_entrez), " / ", length(universe_genes),
  " genes mapped to ENTREZ\n",
  "Modules tested: ", length(module_entrez), "\n",
  "GO ontologies: ", paste(go_ontologies, collapse = ", "), "\n",
  "KEGG: ", if (run_kegg) "yes" else "no", "\n",
  "Cutoffs: pvalue < ", pvalue_cutoff, ", qvalue < ", qvalue_cutoff, "\n"
)

info_path <- file.path(cfg$paths$tables_dir, "05_module_enrichment_run_info.txt")
writeLines(run_info, con = info_path)
message("Wrote run info: ", info_path)

message("Done.")
