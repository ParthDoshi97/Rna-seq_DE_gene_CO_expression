# 08_hub_genes.R
# Turn the kME matrix from 03_wgcna.R into a ranked, trait-anchored,
# DE-cross-referenced hub-gene shortlist. Integrate module preservation
# status from 07_module_preservation.R to identify tumor-specific rewiring.
#
# Flow:
#   - identify which modules are "trait-anchored" (best |cor| & p across all
#     traits, against cfg$hub thresholds)
#   - rank genes within each non-grey module by |kME| in their OWN module column
#   - take top-N per module (optional |kME| floor)
#   - cross-reference DE (log2FoldChange, padj) — non-DE hubs are KEPT;
#     they are the discovery candidates DE alone cannot see
#   - join module preservation status — mark modules not preserved
#     between conditions as biologically significant
#   - sort: priority modules first (trait, not-preserved, or enriched), then
#     by |trait_cor|, then |kME|
#
# Naming stays in "candidate" framing throughout — a hub is a candidate, not
# a driver.

suppressPackageStartupMessages({
  library(yaml)
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

# ---- Inputs ----------------------------------------------------------------

wgcna <- readRDS(cfg$paths$wgcna_rds)
res   <- readRDS(cfg$paths$de_results_rds)
de    <- as.data.frame(res)
de$gene_id <- rownames(res)

moduleColors    <- wgcna$moduleColors
kME             <- wgcna$kME
moduleTraitCor  <- wgcna$moduleTraitCor
moduleTraitPval <- wgcna$moduleTraitPval
datExpr         <- wgcna$datExpr
trait_matrix    <- wgcna$traitMatrix

# kME columns are bare module names (03 strips the "ME" prefix and may drop
# grey). moduleTraitCor rownames are "ME<color>".
stopifnot(length(moduleColors) == ncol(datExpr))

# ---- Config ----------------------------------------------------------------

hub_cfg        <- cfg$hub %||% list()
top_n          <- hub_cfg$top_n          %||% 25
kme_min        <- hub_cfg$kme_min        %||% 0
trait_cor_min  <- hub_cfg$trait_cor_min  %||% 0.3
trait_pval_max <- hub_cfg$trait_pval_max %||% 0.05
enrichment_padj_max <- hub_cfg$enrichment_padj_max %||% 1e-10

# ---- Step A: best trait per module (anchored?) -----------------------------
# For each module row of the module-trait matrix, take the trait with the
# strongest |cor|. Mark the module as trait-anchored if that single best
# pairing clears both thresholds.

mt_modules <- sub("^ME", "", rownames(moduleTraitCor))
trait_anchor <- list()
for (i in seq_along(mt_modules)) {
  m <- mt_modules[i]
  if (m == "grey") next
  cors  <- moduleTraitCor[i, ]
  pvals <- moduleTraitPval[i, ]
  j <- which.max(abs(cors))
  trait_anchor[[m]] <- list(
    best_trait      = colnames(moduleTraitCor)[j],
    best_cor        = unname(cors[j]),
    best_pval       = unname(pvals[j]),
    is_trait_module = isTRUE(abs(cors[j]) >= trait_cor_min &&
                             pvals[j]      <= trait_pval_max)
  )
}

# ---- Step B: rank hubs within each non-grey module -------------------------

modules  <- setdiff(unique(moduleColors), "grey")
hub_rows <- list()
for (m in modules) {
  if (!m %in% colnames(kME)) next               # grey was dropped from kME in 03
  gene_ids <- colnames(datExpr)[moduleColors == m]
  if (length(gene_ids) == 0) next

  kme_vals <- kME[gene_ids, m]
  ord      <- order(abs(kme_vals), decreasing = TRUE)
  ranked   <- gene_ids[ord]
  ranked_k <- kme_vals[ord]

  # Optional absolute |kME| floor — preserves the order of remaining genes.
  keep   <- abs(ranked_k) >= kme_min
  ranked <- ranked[keep]
  ranked_k <- ranked_k[keep]

  take <- min(top_n, length(ranked))
  if (take == 0) next

  hub_rows[[m]] <- data.frame(
    gene_id        = ranked[seq_len(take)],
    module         = m,
    kME            = ranked_k[seq_len(take)],
    rank_in_module = seq_len(take),
    stringsAsFactors = FALSE
  )
}
hub_tbl <- if (length(hub_rows) > 0) {
  do.call(rbind, hub_rows)
} else {
  data.frame(gene_id = character(), module = character(),
             kME = numeric(), rank_in_module = integer())
}
rownames(hub_tbl) <- NULL

message("Ranked hubs: ", nrow(hub_tbl), " across ",
        length(unique(hub_tbl$module)), " module(s)")

# ---- Step C: cross-reference DE -------------------------------------------

hub_tbl$log2FoldChange <- de$log2FoldChange[match(hub_tbl$gene_id, de$gene_id)]
hub_tbl$padj           <- de$padj[match(hub_tbl$gene_id, de$gene_id)]

if ("gene_symbol" %in% colnames(de)) {
  hub_tbl$gene_symbol <- de$gene_symbol[match(hub_tbl$gene_id, de$gene_id)]
} else {
  genes_path <- file.path(dirname(cfg$paths$raw_counts), "genes.tsv")
  if (file.exists(genes_path)) {
    genes <- read.table(genes_path, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
    hub_tbl$gene_symbol <- genes$gene_name[match(hub_tbl$gene_id, genes$gene_id)]
  } else {
    hub_tbl$gene_symbol <- NA_character_
  }
}

padj_cut <- cfg$deseq2$padj_threshold
lfc_cut  <- cfg$deseq2$log2fc_threshold
hub_tbl$is_DE <- !is.na(hub_tbl$padj) & hub_tbl$padj < padj_cut &
                 !is.na(hub_tbl$log2FoldChange) &
                 abs(hub_tbl$log2FoldChange) > lfc_cut

# ---- Step D: attach module-level trait context ----------------------------

ta <- function(m, field, fallback) {
  x <- trait_anchor[[m]]
  if (is.null(x)) fallback else x[[field]]
}
hub_tbl$best_trait <- vapply(hub_tbl$module,
                             ta, "best_trait", NA_character_,
                             FUN.VALUE = character(1))
hub_tbl$trait_cor  <- vapply(hub_tbl$module,
                             ta, "best_cor", NA_real_,
                             FUN.VALUE = numeric(1))
hub_tbl$trait_pval <- vapply(hub_tbl$module,
                             ta, "best_pval", NA_real_,
                             FUN.VALUE = numeric(1))
hub_tbl$is_trait_module <- vapply(hub_tbl$module,
  function(m) isTRUE(trait_anchor[[m]]$is_trait_module),
  logical(1))

# ---- Step D.5: preservation join -------------------------------------------
# Join module preservation status from 07_module_preservation.R.
# Modules NOT preserved between tumor and normal are flagged as biologically
# significant discovery candidates alongside trait-anchored modules.

preservation_file <- cfg$paths$module_preservation_tsv %||%
  file.path(cfg$paths$tables_dir, "module_preservation.tsv")

if (file.exists(preservation_file)) {
  preservation <- read.table(preservation_file, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
  if ("class" %in% colnames(preservation)) {
    pres_map <- setNames(preservation$class, preservation$module)
    hub_tbl$preservation_class <- pres_map[hub_tbl$module]
    hub_tbl$is_not_preserved <- !is.na(hub_tbl$preservation_class) &
      hub_tbl$preservation_class == "not_preserved"
    message("Joined module preservation status: ",
            length(unique(hub_tbl$module[!is.na(hub_tbl$preservation_class)])),
            " modules with preservation data")
  } else {
    hub_tbl$preservation_class <- NA_character_
    hub_tbl$is_not_preserved <- FALSE
    message("Preservation file found but lacks 'class' column; skipping join")
  }
} else {
  hub_tbl$preservation_class <- NA_character_
  hub_tbl$is_not_preserved <- FALSE
  message("No module_preservation.tsv found; priority selection based on ",
          "trait-anchoring alone")
}

# ---- Step D.6: module enrichment join --------------------------------------
# Strong, coherent functional enrichment is a separate discovery signal. This
# catches modules such as greenyellow: preserved in normal tissue, but carrying
# a very strong synaptic/neurotransmitter annotation.

annotation_file <- cfg$paths$module_annotation_summary_csv %||%
  file.path(cfg$paths$tables_dir, "module_annotation_summary.tsv")

hub_tbl$enrichment_best_padj <- NA_real_
hub_tbl$is_enriched_module <- FALSE

if (file.exists(annotation_file)) {
  ann <- read.table(annotation_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, na.strings = c("", "NA"))
  if (all(c("module", "best_padj") %in% colnames(ann))) {
    ann$best_padj <- suppressWarnings(as.numeric(ann$best_padj))
    enrich_map <- setNames(ann$best_padj, ann$module)
    hub_tbl$enrichment_best_padj <- enrich_map[hub_tbl$module]
    hub_tbl$is_enriched_module <- !is.na(hub_tbl$enrichment_best_padj) &
      hub_tbl$enrichment_best_padj <= enrichment_padj_max
    message("Joined module enrichment summary: ",
            length(unique(hub_tbl$module[
              !is.na(hub_tbl$enrichment_best_padj)
            ])),
            " modules with annotation data; ",
            length(unique(hub_tbl$module[hub_tbl$is_enriched_module])),
            " pass enrichment padj <= ", enrichment_padj_max)
  } else {
    message("Module annotation summary found but lacks module/best_padj; ",
            "skipping enrichment priority gate")
  }
} else {
  message("No module_annotation_summary.tsv found; enrichment priority gate ",
          "skipped")
}

# ---- Step E: priority sort -------------------------------------------------
# Define priority modules as trait-anchored OR not-preserved OR strongly
# enriched. Non-DE hubs are deliberately NOT down-ranked.

hub_tbl$is_priority_module <- hub_tbl$is_trait_module |
                              hub_tbl$is_not_preserved |
                              hub_tbl$is_enriched_module

hub_tbl$priority_reason <- vapply(seq_len(nrow(hub_tbl)), function(i) {
  reasons <- character(0)
  if (isTRUE(hub_tbl$is_trait_module[i])) {
    reasons <- c(reasons, "trait")
  }
  if (isTRUE(hub_tbl$is_not_preserved[i])) {
    reasons <- c(reasons, "not_preserved")
  }
  if (isTRUE(hub_tbl$is_enriched_module[i])) {
    reasons <- c(reasons, "enriched")
  }
  if (length(reasons) == 0) "none" else paste(reasons, collapse = ";")
}, character(1))

ord <- order(-as.integer(hub_tbl$is_priority_module),
             -abs(replace(hub_tbl$trait_cor,
                          is.na(hub_tbl$trait_cor), 0)),
             -abs(hub_tbl$kME))
hub_tbl <- hub_tbl[ord, , drop = FALSE]
rownames(hub_tbl) <- NULL

hub_tbl <- hub_tbl[, c("gene_id", "gene_symbol", "module",
                       "kME", "rank_in_module",
                       "log2FoldChange", "padj", "is_DE",
                       "best_trait", "trait_cor", "trait_pval",
                       "is_trait_module",
                       "preservation_class", "is_not_preserved",
                       "enrichment_best_padj", "is_enriched_module",
                       "is_priority_module", "priority_reason")]

# ---- Tables ----------------------------------------------------------------

ranked_path    <- cfg$paths$hub_genes_ranked_csv %||%
  file.path(cfg$paths$tables_dir, "hub_genes_ranked.tsv")
priority_path  <- cfg$paths$hub_genes_priority_csv %||%
  file.path(cfg$paths$tables_dir, "hub_genes_priority_modules.tsv")

write.table(hub_tbl, file = ranked_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(hub_tbl[hub_tbl$is_priority_module, , drop = FALSE],
            file = priority_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote ", nrow(hub_tbl), " ranked hub row(s) -> ", ranked_path)
message("Wrote ", sum(hub_tbl$is_priority_module),
        " priority module hub row(s) ",
        "(trait-anchored, not-preserved, or enriched) -> ", priority_path)

# ---- RDS -------------------------------------------------------------------

hub_rds <- cfg$paths$hub_rds %||%
  file.path(cfg$paths$rds_dir, "hub_genes.rds")
saveRDS(list(
  hub_table    = hub_tbl,
  trait_anchor = trait_anchor,
  thresholds   = list(top_n = top_n, kme_min = kme_min,
                      trait_cor_min = trait_cor_min,
                      trait_pval_max = trait_pval_max,
                      enrichment_padj_max = enrichment_padj_max,
                      padj_threshold = padj_cut,
                      log2fc_threshold = lfc_cut),
  preservation_included = file.exists(preservation_file),
  enrichment_included = file.exists(annotation_file)
), file = hub_rds)

# ---- Figures ---------------------------------------------------------------

hub_fig_dir <- file.path(cfg$paths$figures_dir, "hub_genes")
dir.create(hub_fig_dir, recursive = TRUE, showWarnings = FALSE)

# (a) Barplot of top hubs per priority module — top few modules only.
#     Trait modules are listed first, followed by not-preserved modules and
#     then strongly enriched modules.
anchored_strength <- vapply(trait_anchor,
  function(x) if (isTRUE(x$is_trait_module)) abs(x$best_cor) else NA_real_,
  numeric(1))
anchored_strength <- sort(anchored_strength, decreasing = TRUE, na.last = NA)
trait_modules <- names(anchored_strength)

not_preserved <- setdiff(unique(hub_tbl$module[hub_tbl$is_not_preserved]),
                         trait_modules)
enriched_modules <- unique(hub_tbl$module[hub_tbl$is_enriched_module])
enriched_padj <- vapply(enriched_modules, function(m) {
  vals <- hub_tbl$enrichment_best_padj[hub_tbl$module == m]
  min(vals, na.rm = TRUE)
}, numeric(1))
enriched_modules <- enriched_modules[order(enriched_padj)]
enriched_modules <- setdiff(enriched_modules, c(trait_modules, not_preserved))

max_plot_modules <- hub_cfg$plot_max_priority_modules %||% 12
top_modules_for_bar <- unique(c(trait_modules, not_preserved, enriched_modules))
top_modules_for_bar <- head(top_modules_for_bar, max_plot_modules)

for (m in top_modules_for_bar) {
  sub <- hub_tbl[hub_tbl$module == m, , drop = FALSE]
  if (nrow(sub) == 0) next
  sub <- sub[order(abs(sub$kME), decreasing = TRUE), , drop = FALSE]
  labels <- ifelse(is.na(sub$gene_symbol) | sub$gene_symbol == "",
                   sub$gene_id, sub$gene_symbol)
  best_cor <- trait_anchor[[m]]$best_cor
  best_tr  <- trait_anchor[[m]]$best_trait

  png(file.path(hub_fig_dir, paste0("hubs_", m, ".png")),
      width = 8, height = max(4, 0.25 * nrow(sub) + 2),
      units = "in", res = 150)
  par(mar = c(4, 9, 4, 2))
  barplot(rev(abs(sub$kME)), names.arg = rev(labels),
          horiz = TRUE, las = 1,
          col = m, border = NA,
          xlab = "|kME| (module membership)",
          main = paste0("Hub candidates - module ", m,
                        " (", best_tr, ", r=", round(best_cor, 2), ")"))
  dev.off()
}

# (b) GS-vs-MM scatter per priority module — the classic plot. Hubs sit in
#     the top-right: high module membership AND high gene-trait significance.
gs_dir <- file.path(hub_fig_dir, "gs_vs_mm")
dir.create(gs_dir, recursive = TRUE, showWarnings = FALSE)
for (m in top_modules_for_bar) {
  trait_name <- trait_anchor[[m]]$best_trait
  if (is.null(trait_name) || !trait_name %in% colnames(trait_matrix)) next
  in_module <- moduleColors == m
  gs <- abs(drop(cor(datExpr, trait_matrix[, trait_name], use = "p")))
  mm <- abs(kME[, m])

  png(file.path(gs_dir,
                paste0("gs_vs_mm_", m, "_", make.names(trait_name), ".png")),
      width = 7, height = 7, units = "in", res = 150)
  par(mar = c(5, 5, 4, 2))
  plot(mm[in_module], gs[in_module],
       xlab = paste0("|kME| in module ", m),
       ylab = paste0("|cor(gene, ", trait_name, ")|"),
       main = paste0("GS vs MM - ", m, " / ", trait_name),
       pch = 20, col = m)
  abline(lm(gs[in_module] ~ mm[in_module]), col = "black", lty = 2)
  dev.off()
}

message("Done.")
