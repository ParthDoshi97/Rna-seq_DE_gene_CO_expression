#!/usr/bin/env Rscript
# 03_wgcna.R
# WGCNA on the variance-stabilised expression matrix:
# pick soft-threshold → adjacency / TOM → dynamic-tree-cut modules
# → merge close modules → module–trait correlation → kME hub candidates.
# (Branch 5 adds the scale-free-topology fit plot, sample-clustering tree,
#  and the per-module hub-gene deliverable.)

suppressPackageStartupMessages({
  library(yaml)
  library(SummarizedExperiment)
  library(DESeq2)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
cfg <- yaml::read_yaml(config_path)

set.seed(cfg$runtime$seed)
WGCNA::allowWGCNAThreads(nThreads = cfg$runtime$threads)
dir.create(cfg$paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$rds_dir,     recursive = TRUE, showWarnings = FALSE)

# ---- Input: variance-stabilised matrix from script 02 -----------------------

vsd <- readRDS(cfg$paths$vsd_rds)
samples <- as.data.frame(colData(vsd))
samples$condition <- factor(
  samples$condition,
  levels = c(cfg$dataset$reference_level, cfg$dataset$treatment_level)
)

expr_full <- assay(vsd)
vars <- matrixStats::rowVars(expr_full)
keep <- vars > quantile(vars, cfg$wgcna$variance_quantile)
message("WGCNA gene filter: kept ", sum(keep), " / ", length(keep),
        " genes above variance quantile ", cfg$wgcna$variance_quantile)
datExpr <- t(expr_full[keep, , drop = FALSE])  # samples × genes

# ---- Soft-threshold power ---------------------------------------------------

sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = 1:20,
  networkType = cfg$wgcna$network_type,
  verbose = 0
)
soft_power <- cfg$wgcna$soft_power
message("Soft power used: ", soft_power,
        "  (pickSoftThreshold suggested: ",
        ifelse(is.na(sft$powerEstimate), "NA", sft$powerEstimate), ")")

# ---- Adjacency + TOM + modules ----------------------------------------------

adj <- WGCNA::adjacency(datExpr, power = soft_power,
                        type = cfg$wgcna$network_type)
TOM <- WGCNA::TOMsimilarity(adj, TOMType = cfg$wgcna$tom_type)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

dynamicMods <- cutreeDynamic(
  dendro = geneTree, distM = dissTOM,
  deepSplit = cfg$wgcna$deep_split,
  pamRespectsDendro = cfg$wgcna$pam_respects_dendro,
  minClusterSize = cfg$wgcna$min_module_size
)
dynamicColors <- WGCNA::labels2colors(dynamicMods)

# Merge modules whose eigengenes are very similar
merged <- WGCNA::mergeCloseModules(
  datExpr, dynamicColors,
  cutHeight = cfg$wgcna$merge_cut_height, verbose = 0
)
moduleColors <- merged$colors
MEs <- merged$newMEs

message("Modules after merge: ",
        length(unique(moduleColors)),
        "  (sizes: ",
        paste(sort(table(moduleColors), decreasing = TRUE), collapse = ", "),
        ")")

# ---- Module–trait correlation -----------------------------------------------

design <- model.matrix(~ condition - 1, data = samples)
colnames(design) <- levels(samples$condition)

moduleTraitCor   <- cor(MEs, design, use = "p")
moduleTraitPval  <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(datExpr))

# ---- Hub gene candidacy by kME ----------------------------------------------

kME <- WGCNA::signedKME(datExpr, MEs)
colnames(kME) <- sub("^kME", "", colnames(kME))

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
    moduleColors = moduleColors,
    MEs = MEs,
    moduleTraitCor = moduleTraitCor,
    moduleTraitPval = moduleTraitPval,
    kME = kME,
    geneTree = geneTree,
    soft_power = soft_power,
    sft = sft,
    datExpr = datExpr
  ),
  file = cfg$paths$wgcna_rds
)

# ---- Core figures (others added in branch 5) --------------------------------

png(file.path(cfg$paths$figures_dir, "wgcna_dendrogram.png"),
    width = 12, height = 8, units = "in", res = 150)
WGCNA::plotDendroAndColors(
  geneTree, moduleColors,
  groupLabels = "Module", dendroLabels = FALSE,
  hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = paste0("Gene clustering with module colours — ",
                cfg$project$cancer_short)
)
dev.off()

png(file.path(cfg$paths$figures_dir, "wgcna_module_trait.png"),
    width = 8, height = 10, units = "in", res = 150)
par(mar = c(6, 9, 4, 2))
textMat <- paste0(signif(moduleTraitCor, 2),
                  "\n(", signif(moduleTraitPval, 1), ")")
dim(textMat) <- dim(moduleTraitCor)
WGCNA::labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(design), yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE, colors = WGCNA::blueWhiteRed(50),
  textMatrix = textMat, setStdMargins = FALSE,
  cex.text = 0.6, zlim = c(-1, 1),
  main = paste0("Module–trait relationships — ", cfg$project$cancer_short)
)
dev.off()
