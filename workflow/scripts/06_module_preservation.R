# 07_module_preservation.R
# Module preservation between the tumor (reference) and normal (test)
# networks. A tumor module that is NOT preserved in normal samples indicates
# a co-expression program specific to the tumor state — i.e. network rewiring
# that DE alone cannot see.
#
# This is the only stage that re-includes the normal arm: 03_wgcna.R was run
# tumor-only by design so modules reflect tumor-internal variation, but
# preservation requires both arms to make the comparison.
#
# Flow:
#   - load tumor WGCNA results and module assignments (from 03)
#   - load all (tumor + normal) count data and expression
#   - run modulePreservation(testSet, referenceSet, ...)
#   - output: module_preservation.tsv with Zsummary, medianRank, class
#   - class: "preserved" (Z > 5), "moderate" (2 < Z < 5), "not_preserved" (Z < 2)

suppressPackageStartupMessages({
  library(yaml)
  library(WGCNA)
  library(dplyr)
  library(SummarizedExperiment)
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
allowWGCNAThreads(nThreads = cfg$runtime$threads)
dir.create(cfg$paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$rds_dir,     recursive = TRUE, showWarnings = FALSE)

# ---- Inputs ----------------------------------------------------------------

wgcna_tumor <- readRDS(cfg$paths$wgcna_rds)
module_assignments_path <- cfg$paths$module_assignments_csv %||%
  file.path(cfg$paths$tables_dir, "module_assignments.tsv")
module_assignments <- read.table(module_assignments_path,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# Load full dataset (tumor + normal)
counts <- read.table(cfg$paths$raw_counts, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, row.names = 1)
samples <- read.table(cfg$paths$sample_metadata, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, row.names = 1)

# Config
pres_cfg <- cfg$preservation %||% list()
n_permutations <- pres_cfg$n_permutations %||% 200
min_normal <- pres_cfg$min_normal %||% 20
max_module_size <- pres_cfg$max_module_size %||% 2000

condition_col <- cfg$dataset$condition_column
ref_level <- cfg$dataset$reference_level
treat_level <- cfg$dataset$treatment_level

# ---- Prepare expression data ------------------------------------------------

ref_idx <- rownames(samples)[samples[[condition_col]] == ref_level]
test_idx <- rownames(samples)[samples[[condition_col]] == treat_level]

message("Reference (normal): ", length(ref_idx), " samples")
message("Test (tumor): ", length(test_idx), " samples")

if (length(ref_idx) < min_normal) {
  stop("Fewer than min_normal reference samples: ", length(ref_idx), " < ", min_normal)
}

# Filter counts to variance-selected genes (match WGCNA input)
vsd_rds <- cfg$paths$vsd_rds
vsd <- readRDS(vsd_rds)
expr_full <- t(assay(vsd))

# Split by condition: rows = samples, cols = genes
expr_ref <- expr_full[rownames(expr_full) %in% ref_idx, ]
expr_test <- expr_full[rownames(expr_full) %in% test_idx, ]

message("Expression (reference): ", nrow(expr_ref), " samples x ",
        ncol(expr_ref), " genes")
message("Expression (test): ", nrow(expr_test), " samples x ",
        ncol(expr_test), " genes")

# ---- Filter expression data for quality ------------------------------------
# modulePreservation requires genes/samples without zero variance or excessive NAs

message("Filtering expression data for quality...")

# Check test set
gsg_test <- WGCNA::goodSamplesGenes(expr_test, verbose = 1)
if (!gsg_test$allOK) {
  message("  Test set: removing ", sum(!gsg_test$goodGenes), " bad genes, ",
          sum(!gsg_test$goodSamples), " bad samples")
  expr_test <- expr_test[gsg_test$goodSamples, gsg_test$goodGenes]
}

# Check reference set
gsg_ref <- WGCNA::goodSamplesGenes(expr_ref, verbose = 1)
if (!gsg_ref$allOK) {
  message("  Reference set: removing ", sum(!gsg_ref$goodGenes), " bad genes, ",
          sum(!gsg_ref$goodSamples), " bad samples")
  expr_ref <- expr_ref[gsg_ref$goodSamples, gsg_ref$goodGenes]
}

# Ensure both networks have identical gene sets for modulePreservation
# (find intersection of good genes from both networks)
common_genes <- intersect(colnames(expr_test), colnames(expr_ref))
if (length(common_genes) < ncol(expr_test)) {
  message("  Filtering to common genes: ", length(common_genes), " / ",
          ncol(expr_test), " from test set")
  expr_test <- expr_test[, common_genes]
}
if (length(common_genes) < ncol(expr_ref)) {
  message("  Filtering to common genes: ", length(common_genes), " / ",
          ncol(expr_ref), " from reference set")
  expr_ref <- expr_ref[, common_genes]
}

# Preservation is only defined for genes that were assigned to modules by
# 03_wgcna.R. The VST object contains every filtered expression gene, so keep
# the network input genes here instead of passing thousands of non-module genes
# and relying on WGCNA to remove missing colors later.
network_genes <- intersect(colnames(wgcna_tumor$datExpr), common_genes)
if (length(network_genes) == 0) {
  stop("No overlap between WGCNA module genes and preservation expression data.")
}
message("  Filtering to WGCNA module genes: ", length(network_genes), " / ",
        length(common_genes), " common genes")
expr_test <- expr_test[, network_genes, drop = FALSE]
expr_ref  <- expr_ref[,  network_genes, drop = FALSE]

message("After filtering:")
message("  Reference: ", nrow(expr_ref), " samples x ", ncol(expr_ref), " genes")
message("  Test: ", nrow(expr_test), " samples x ", ncol(expr_test), " genes")

# ---- Prepare module assignments for preservation ---------------------------

# Build a vector: names = gene_id, values = numeric module label.
# modulePreservation() wants numbers, but our downstream tables use WGCNA color
# names. Keep a label -> color map so results can be joined back to hub genes.
module_index <- match(colnames(expr_test), module_assignments$gene_id)
if (anyNA(module_index)) {
  stop("Missing module assignments for ",
       sum(is.na(module_index)), " preservation gene(s).")
}
module_levels <- c("grey", setdiff(unique(module_assignments$module), "grey"))
module_label_map <- setNames(module_levels, as.character(seq_along(module_levels) - 1L))
module_label_map["0.1"] <- "all_modules"

module_vector <- as.numeric(factor(
  module_assignments$module[module_index],
  levels = module_levels
)) - 1L
names(module_vector) <- colnames(expr_test)

# Log module sizes
module_colors <- unique(sort(module_assignments$module))
message("\nModule sizes:")
for (m in module_colors) {
  size <- sum(module_assignments$module == m)
  if (size > max_module_size) {
    message("  ", m, ": ", size, " (WARNING: exceeds max_module_size ",
            max_module_size, ")")
  } else {
    message("  ", m, ": ", size)
  }
}

# ---- Run preservation analysis -----------------------------------------------

message("\nRunning modulePreservation with ", n_permutations, " permutations...")

# modulePreservation requires: multiData = list of networks with data
# Network 1 = test (tumor), Network 2 = reference (normal)
preservation <- modulePreservation(
  multiData = list(
    list(data = expr_test),  # Network 1: test (tumor)
    list(data = expr_ref)    # Network 2: reference (normal)
  ),
  multiColor = list(module_vector, module_vector),
  referenceNetworks = 2,
  testNetworks = 1,
  nPermutations = n_permutations,
  maxModuleSize = max_module_size,
  verbose = 3
)

# ---- Extract and format results --------------------------------------------

# WGCNA modulePreservation() output structure can vary by WGCNA version.
# Some versions use:
#   preservation$preservation$Z$ref.2$inColumnsAlsoPresentIn.1
# Other versions use:
#   preservation$preservation$Z$ref.Set_2$inColumnsAlsoPresentIn.Set_1

message("Preservation result top-level names: ",
        paste(names(preservation), collapse = ", "))

pres_inner <- preservation$preservation
if (is.null(pres_inner)) {
  pres_inner <- preservation
}

message("Preservation inner names: ",
        paste(names(pres_inner), collapse = ", "))

message("Available Z reference paths: ",
        paste(names(pres_inner$Z), collapse = ", "))

message("Available observed reference paths: ",
        paste(names(pres_inner$observed), collapse = ", "))

# Robustly find reference and test result paths
z_ref_name <- grep("^ref", names(pres_inner$Z), value = TRUE)[1]
obs_ref_name <- grep("^ref", names(pres_inner$observed), value = TRUE)[1]

if (is.na(z_ref_name) || is.na(obs_ref_name)) {
  stop("Could not locate reference paths in preservation result.")
}

message("Available Z test paths under ", z_ref_name, ": ",
        paste(names(pres_inner$Z[[z_ref_name]]), collapse = ", "))

message("Available observed test paths under ", obs_ref_name, ": ",
        paste(names(pres_inner$observed[[obs_ref_name]]), collapse = ", "))

z_test_name <- grep("inColumnsAlsoPresentIn",
                    names(pres_inner$Z[[z_ref_name]]),
                    value = TRUE)[1]

obs_test_name <- grep("inColumnsAlsoPresentIn",
                      names(pres_inner$observed[[obs_ref_name]]),
                      value = TRUE)[1]

if (is.na(z_test_name) || is.na(obs_test_name)) {
  stop("Could not locate test-network paths in preservation result.")
}

message("Using Z path: Z$", z_ref_name, "$", z_test_name)
message("Using observed path: observed$", obs_ref_name, "$", obs_test_name)

z_stats <- pres_inner$Z[[z_ref_name]][[z_test_name]]
obs_stats <- pres_inner$observed[[obs_ref_name]][[obs_test_name]]

if (is.null(z_stats)) {
  stop("Could not locate Z statistics in preservation result.")
}

message("Z stats columns: ", paste(colnames(z_stats), collapse = ", "))
message("Z stats rows/modules: ", paste(rownames(z_stats), collapse = ", "))

# Handle possible column-name differences safely
zsummary_col <- grep("^Zsummary", colnames(z_stats), value = TRUE)[1]
medianrank_col <- grep("^medianRank", colnames(z_stats), value = TRUE)[1]

if (is.na(zsummary_col)) {
  stop("Could not find Zsummary column. Available columns: ",
       paste(colnames(z_stats), collapse = ", "))
}

if (is.na(medianrank_col)) {
  warning("Could not find medianRank column. Setting medianRank to NA.")
}

raw_modules <- rownames(z_stats)
module_names <- unname(module_label_map[raw_modules])
module_names[is.na(module_names)] <- raw_modules[is.na(module_names)]

pres_df <- data.frame(
  module_label = raw_modules,
  module = module_names,
  moduleSize = NA_integer_,
  Zsummary = z_stats[[zsummary_col]],
  medianRank = if (!is.na(medianrank_col)) {
    z_stats[[medianrank_col]]
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)

# Classify preservation
pres_df$class <- ifelse(
  is.na(pres_df$Zsummary), "unknown",
  ifelse(
    pres_df$Zsummary > 5, "preserved",
    ifelse(pres_df$Zsummary > 2, "moderate", "not_preserved")
  )
)

# Add module sizes from the original color assignment table. The synthetic
# "all_modules" row is useful as a WGCNA diagnostic but is not a real module.
module_sizes <- table(module_assignments$module)
pres_df$moduleSize <- as.integer(module_sizes[pres_df$module])
pres_df$moduleSize[pres_df$module == "all_modules"] <- sum(module_sizes)

# Sort by Zsummary descending
pres_df <- pres_df[order(-pres_df$Zsummary, na.last = TRUE), ]
rownames(pres_df) <- NULL

message("\nModule preservation summary:")
print(pres_df)

# ---- Outputs ---------------------------------------------------------------

preservation_path <- cfg$paths$module_preservation_tsv %||%
  file.path(cfg$paths$tables_dir, "module_preservation.tsv")
write.table(pres_df, file = preservation_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote preservation table: ", preservation_path)

preservation_rds <- cfg$paths$preservation_rds %||%
  file.path(cfg$paths$rds_dir, "preservation.rds")
saveRDS(list(
  preservation = preservation,
  table = pres_df,
  module_label_map = module_label_map,
  n_permutations = n_permutations
), file = preservation_rds)
message("Wrote preservation RDS: ", preservation_rds)

# ---- Figures ---------------------------------------------------------------

fig_dir <- file.path(cfg$paths$figures_dir, "07_module_preservation")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Plot Zsummary per module
tryCatch({
  pres_plot <- pres_df[!pres_df$module %in% c("grey", "all_modules"), ]
  if (nrow(pres_plot) == 0) {
    stop("no non-grey modules to plot")
  }
  pres_plot <- pres_plot[order(pres_plot$Zsummary), ]
  pres_plot$module <- factor(pres_plot$module,
                             levels = pres_plot$module)

  png(file.path(fig_dir, "preservation_zsummary.png"),
      width = 8, height = max(4, 0.3 * nrow(pres_plot) + 2),
      units = "in", res = 150)
  par(mar = c(5, 8, 4, 2))
  col_map <- setNames(
    c("darkred", "orange", "darkgreen"),
    c("not_preserved", "moderate", "preserved")
  )
  cols <- col_map[pres_plot$class]
  barplot(pres_plot$Zsummary, names.arg = pres_plot$module,
          horiz = TRUE, las = 1, col = cols, border = NA,
          xlab = "Z-summary (preservation)")
  abline(v = c(2, 5), lty = 2, col = "grey")
  legend("bottomright", c("Not preserved", "Moderate", "Preserved"),
         fill = c("darkred", "orange", "darkgreen"), bty = "n")
  dev.off()
  message("Wrote figure: ", fig_dir, "/preservation_zsummary.png")
}, error = function(e) {
  warning("Zsummary plot failed: ", e$message)
})

# Plot median rank per module
tryCatch({
  pres_plot <- pres_df[!pres_df$module %in% c("grey", "all_modules") &
                         !is.na(pres_df$medianRank), ]
  if (nrow(pres_plot) == 0) {
    message("Median rank plot skipped: no finite medianRank values returned")
  } else {
    pres_plot <- pres_plot[order(pres_plot$medianRank), ]
    pres_plot$module <- factor(pres_plot$module,
                               levels = pres_plot$module)

    png(file.path(fig_dir, "preservation_medianrank.png"),
        width = 8, height = max(4, 0.3 * nrow(pres_plot) + 2),
        units = "in", res = 150)
    par(mar = c(5, 8, 4, 2))
    barplot(pres_plot$medianRank, names.arg = pres_plot$module,
            horiz = TRUE, las = 1, col = pres_plot$module, border = NA,
            xlab = "Median rank (preservation)")
    dev.off()
    message("Wrote figure: ", fig_dir, "/preservation_medianrank.png")
  }
}, error = function(e) {
  warning("Median rank plot failed: ", e$message)
})

# Write run info
run_info <- paste0(
  "Module Preservation Run Info\n",
  "=============================\n",
  "Timestamp: ", Sys.time(), "\n",
  "Script: 07_module_preservation.R\n",
  "WGCNA version: ", packageVersion("WGCNA"), "\n",
  "\n",
  "Reference (normal): ", length(ref_idx), " samples\n",
  "Test (tumor): ", length(test_idx), " samples\n",
  "Permutations: ", n_permutations, "\n",
  "\n",
  "Classification:\n",
  "  Z > 5: preserved\n",
  "  2 < Z < 5: moderate\n",
  "  Z < 2: not_preserved\n"
)

info_path <- file.path(cfg$paths$tables_dir, "07_module_preservation_run_info.txt")
writeLines(run_info, con = info_path)
message("Wrote run info: ", info_path)

message("Done.")
