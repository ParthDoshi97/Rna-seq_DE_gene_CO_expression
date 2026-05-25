# 03_wgcna.R
# Tumor-only WGCNA on the variance-stabilised expression matrix.
# Flow:
#   sample QC (goodSamplesGenes + sample-tree outlier cut)
#   -> pick soft-threshold (with scale-free R^2 check at the chosen power)
#   -> adjacency / TOM -> dynamic-tree-cut modules -> merge close modules
#   -> module-trait correlation against clinical traits
#      (ER/PR/HER2/stage/age, not the tumor-vs-normal axis)
#   -> kME hub-gene candidates (grey module excluded).

suppressPackageStartupMessages({
  library(yaml)
  library(SummarizedExperiment)
  library(DESeq2)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
cfg <- yaml::read_yaml(config_path)

set.seed(cfg$runtime$seed)
WGCNA::allowWGCNAThreads(nThreads = cfg$runtime$threads)
dir.create(cfg$paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$rds_dir,     recursive = TRUE, showWarnings = FALSE)

# ---- Consistency guard ------------------------------------------------------
# Signed networks need signed TOM. Catching this here beats discovering it
# later as a one-giant-module surprise.
if (cfg$wgcna$network_type %in% c("signed", "signed hybrid") &&
    !grepl("^signed", cfg$wgcna$tom_type)) {
  stop("network_type='", cfg$wgcna$network_type,
       "' requires tom_type to start with 'signed' (got '",
       cfg$wgcna$tom_type, "')")
}

# ---- Input: variance-stabilised matrix from script 02 -----------------------

vsd <- readRDS(cfg$paths$vsd_rds)
samples <- as.data.frame(colData(vsd))
expr_full <- assay(vsd)

# ---- Tumor-only filter ------------------------------------------------------
# With both arms in the matrix the first principal component IS the contrast,
# and modules collapse to "cancer vs not". Filter to the treatment level
# BEFORE the variance cut so the kept genes reflect tumor-internal variation.

if (isTRUE(cfg$wgcna$tumor_only)) {
  tumor_keep <- samples$condition == cfg$dataset$treatment_level
  message("Tumor-only filter: keeping ", sum(tumor_keep), " / ",
          length(tumor_keep),
          " samples (condition == '", cfg$dataset$treatment_level, "')")
  expr_full <- expr_full[, tumor_keep, drop = FALSE]
  samples   <- samples[tumor_keep, , drop = FALSE]
}

# ---- Variance filter on the post-tumor-only matrix --------------------------
# Top-N by variance rather than a quantile: keeps the gene count deterministic
# regardless of how many genes survived prefilter. Too many genes leaks the
# low-variance tail into the grey/unassigned module.

vars  <- matrixStats::rowVars(expr_full)
n_top <- min(cfg$wgcna$n_top_genes %||% 5000, length(vars))
keep  <- rank(-vars, ties.method = "first") <= n_top
message("WGCNA gene filter: kept top ", sum(keep), " / ", length(keep),
        " genes by variance")
datExpr <- t(expr_full[keep, , drop = FALSE])     # samples x genes

# ---- WGCNA sample/gene QC ---------------------------------------------------

gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  message("goodSamplesGenes dropped ",
          sum(!gsg$goodGenes),  " gene(s), ",
          sum(!gsg$goodSamples), " sample(s)")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  samples <- samples[gsg$goodSamples, , drop = FALSE]
}

# ---- Sample outlier removal via hierarchical clustering ---------------------

sampleTree <- hclust(dist(datExpr), method = "average")
sample_cut <- cfg$wgcna$sample_cut_height
if (!is.null(sample_cut) && is.numeric(sample_cut)) {
  clust    <- WGCNA::cutreeStatic(sampleTree, cutHeight = sample_cut,
                                  minSize = 10)
  keepSamp <- clust == 1
  n_drop   <- sum(!keepSamp)
  if (n_drop > 0) {
    message("Sample outlier cut (cutHeight=", sample_cut, "): dropping ",
            n_drop, " sample(s)")
    datExpr <- datExpr[keepSamp, , drop = FALSE]
    samples <- samples[keepSamp, , drop = FALSE]
  } else {
    message("Sample outlier cut: no outliers above cutHeight=", sample_cut)
  }
} else {
  message("Sample outlier cut: skipped (wgcna.sample_cut_height is null). ",
          "Inspect wgcna_sample_tree.png and set a cutoff.")
}

# Sample-tree figure (always emit -- it's the visual companion to the cut).
png(file.path(cfg$paths$figures_dir, "wgcna_sample_tree.png"),
    width = 12, height = 6, units = "in", res = 150)
par(mar = c(2, 4, 4, 2))
plot(sampleTree, labels = FALSE,
     main = paste0("Sample clustering -- ", cfg$project$cancer_short),
     sub = "", xlab = "")
if (!is.null(sample_cut) && is.numeric(sample_cut)) {
  abline(h = sample_cut, col = "red", lty = 2)
}
dev.off()

# ---- Soft-threshold power + scale-free R^2 check ----------------------------

sft <- WGCNA::pickSoftThreshold(
  datExpr, powerVector = 1:20,
  networkType = cfg$wgcna$network_type, verbose = 0
)

soft_power <- cfg$wgcna$soft_power
if (is.null(soft_power)) {
  if (!is.na(sft$powerEstimate)) {
    soft_power <- sft$powerEstimate
    message("soft_power not in config; using pickSoftThreshold estimate: ",
            soft_power)
  } else {
    stop("No soft_power in config and pickSoftThreshold returned NA. ",
         "Pick a power manually after inspecting wgcna_soft_threshold.png.")
  }
}

fit         <- sft$fitIndices
r2_at_power <- fit$SFT.R.sq[fit$Power == soft_power]
r2_min      <- cfg$wgcna$scale_free_r2_min %||% 0.80
message("Soft power used: ", soft_power,
        "  (pickSoftThreshold suggested: ",
        ifelse(is.na(sft$powerEstimate), "NA", sft$powerEstimate), ")",
        "  scale-free R^2 = ",
        ifelse(length(r2_at_power) == 0, "NA", round(r2_at_power, 3)))
if (length(r2_at_power) == 0 || r2_at_power < r2_min) {
  warning("Scale-free R^2 at power ", soft_power, " (",
          ifelse(length(r2_at_power) == 0, "NA", round(r2_at_power, 3)),
          ") < threshold ", r2_min,
          " -- network may be noise-dominated. Re-tune wgcna.soft_power.")
}

# Soft-threshold diagnostic plot.
png(file.path(cfg$paths$figures_dir, "wgcna_soft_threshold.png"),
    width = 10, height = 5, units = "in", res = 150)
par(mfrow = c(1, 2))
plot(fit$Power, -sign(fit$slope) * fit$SFT.R.sq,
     xlab = "Soft power", ylab = "Scale-free fit R^2",
     type = "n", main = "Scale independence")
text(fit$Power, -sign(fit$slope) * fit$SFT.R.sq,
     labels = fit$Power, col = "red")
abline(h = r2_min, col = "blue", lty = 2)
plot(fit$Power, fit$mean.k.,
     xlab = "Soft power", ylab = "Mean connectivity",
     type = "n", main = "Mean connectivity")
text(fit$Power, fit$mean.k., labels = fit$Power, col = "red")
dev.off()

# ---- Adjacency + TOM + modules ----------------------------------------------

adj <- WGCNA::adjacency(datExpr, power = soft_power,
                        type = cfg$wgcna$network_type)
TOM <- WGCNA::TOMsimilarity(adj, TOMType = cfg$wgcna$tom_type)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

dynamicMods <- cutreeDynamic(
  dendro = geneTree, distM = dissTOM,
  deepSplit = cfg$wgcna$deep_split,
  pamStage = TRUE,   # PAM reassigns borderline genes off the grey pile to their nearest module
  pamRespectsDendro = cfg$wgcna$pam_respects_dendro,
  minClusterSize = cfg$wgcna$min_module_size
)
dynamicColors <- WGCNA::labels2colors(dynamicMods)

# Merge modules whose eigengenes are very similar.
merged <- WGCNA::mergeCloseModules(
  datExpr, dynamicColors,
  cutHeight = cfg$wgcna$merge_cut_height, verbose = 0
)
moduleColors <- merged$colors
MEs <- merged$newMEs

message("Modules after merge: ", length(unique(moduleColors)),
        "  (sizes: ",
        paste(sort(table(moduleColors), decreasing = TRUE), collapse = ", "),
        ")")

# ---- Module-trait correlation using clinical traits -------------------------
# Once tumor_only is on, the condition column is constant and useless as a
# trait. 01_load_and_qc.R has copied a configurable set of clinical columns
# (ER/PR/HER2/stage/...) from recount3 colData into samples.tsv. Pick which
# ones to correlate via cfg$traits$module_trait_columns.

build_trait_matrix <- function(samples, cols) {
  na_tokens <- c("", "NA", "[Not Available]", "[Unknown]",
                 "[Not Evaluated]", "[Not Applicable]",
                 "Indeterminate", "Equivocal")
  out_list <- list()
  for (col in cols) {
    if (!col %in% colnames(samples)) {
      message("Trait '", col, "' not in samples; skipping")
      next
    }
    v <- samples[[col]]
    if (is.numeric(v) || is.logical(v)) {
      out_list[[col]] <- as.numeric(v)
    } else {
      vc <- as.character(v)
      vc[vc %in% na_tokens] <- NA
      lev <- sort(unique(vc[!is.na(vc)]))
      if (length(lev) < 2) {
        message("Trait '", col, "' has <2 informative levels; skipping")
        next
      }
      if (length(lev) == 2) {
        # 0/1 binary against the alphabetically-first level.
        out_list[[col]] <- as.integer(factor(vc, levels = lev)) - 1L
      } else {
        # One-hot for multi-level categoricals (e.g. stage).
        for (l in lev) {
          out_list[[paste0(col, "_", l)]] <- as.integer(vc == l)
        }
      }
    }
  }
  if (length(out_list) == 0) return(NULL)
  out <- do.call(cbind, lapply(out_list, function(x) x))
  rownames(out) <- samples$sample_id
  out
}

trait_cols   <- cfg$traits$module_trait_columns %||% character(0)
trait_matrix <- build_trait_matrix(samples, trait_cols)
if (is.null(trait_matrix) || ncol(trait_matrix) == 0) {
  warning("No usable trait columns -- falling back to lib_size-only design.")
  fallback <- if ("lib_size" %in% colnames(samples)) {
    as.numeric(samples$lib_size)
  } else {
    colSums(datExpr)
  }
  trait_matrix <- matrix(fallback, ncol = 1,
                         dimnames = list(samples$sample_id, "lib_size"))
}
# Hard guard: trait correlations are only meaningful if rows align by sample ID.
# Assigning by position (the prior bug) would silently miscorrelate every trait
# if anything ever reordered `samples` relative to `datExpr` upstream.
stopifnot(identical(rownames(trait_matrix), rownames(datExpr)))
message("Module-trait matrix: ", ncol(trait_matrix), " trait(s) (",
        paste(colnames(trait_matrix), collapse = ", "), ")")

# MEgrey is the unassigned bin, not a co-expression module — including it in the
# trait correlation produces a misleading row in the heatmap. Drop it here but
# keep `MEs` intact for downstream kME computation.
MEs_for_trait   <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
moduleTraitCor  <- cor(MEs_for_trait, trait_matrix, use = "p")
moduleTraitPval <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(datExpr))

# ---- Hub gene candidacy by kME (grey excluded) ------------------------------

kME <- WGCNA::signedKME(datExpr, MEs)
colnames(kME) <- sub("^kME", "", colnames(kME))
if (isTRUE(cfg$wgcna$exclude_grey_in_hub %||% TRUE) &&
    "grey" %in% colnames(kME)) {
  kME <- kME[, colnames(kME) != "grey", drop = FALSE]
  message("kME: dropped 'grey' (unassigned bin) from hub-gene candidate pool")
}

# ---- Deliverables -----------------------------------------------------------

module_tbl <- data.frame(
  gene_id = colnames(datExpr),
  module  = moduleColors,
  stringsAsFactors = FALSE
)
write.table(module_tbl, file = cfg$paths$module_assignments_csv,
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(
  list(
    moduleColors    = moduleColors,
    MEs             = MEs,
    moduleTraitCor  = moduleTraitCor,
    moduleTraitPval = moduleTraitPval,
    traitMatrix     = trait_matrix,
    kME             = kME,
    geneTree        = geneTree,
    sampleTree      = sampleTree,
    soft_power      = soft_power,
    sft             = sft,
    datExpr         = datExpr,
    samples         = samples
  ),
  file = cfg$paths$wgcna_rds
)

# ---- Core figures -----------------------------------------------------------

png(file.path(cfg$paths$figures_dir, "wgcna_dendrogram.png"),
    width = 12, height = 8, units = "in", res = 150)
WGCNA::plotDendroAndColors(
  geneTree, moduleColors,
  groupLabels = "Module", dendroLabels = FALSE,
  hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = paste0("Gene clustering with module colours -- ",
                cfg$project$cancer_short)
)
dev.off()

# Module-trait heatmap — width scales with trait count.
trait_w <- max(6, 0.6 * ncol(moduleTraitCor) + 3)
png(file.path(cfg$paths$figures_dir, "wgcna_module_trait.png"),
    width = trait_w, height = 10, units = "in", res = 150)
par(mar = c(8, 9, 4, 2))
textMat <- paste0(signif(moduleTraitCor, 2),
                  "\n(", signif(moduleTraitPval, 1), ")")
dim(textMat) <- dim(moduleTraitCor)
WGCNA::labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(moduleTraitCor),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE, colors = WGCNA::blueWhiteRed(50),
  textMatrix = textMat, setStdMargins = FALSE,
  cex.text = 0.6, zlim = c(-1, 1),
  main = paste0("Module-trait relationships -- ", cfg$project$cancer_short)
)
dev.off()
