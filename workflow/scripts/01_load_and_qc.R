# 01_load_and_qc.R
# Pull TCGA-BRCA counts via recount3 (or load a local counts table),
# build the sample table, apply biotype + low-count gene filtering and
# sample-level QC, and emit counts.tsv + samples.tsv + genes.tsv that
# downstream scripts consume.

# ---- Install missing packages (CRAN + Bioconductor) -------------------------
# Fallback for environments without environment.yml / the container.
# The pinned environment (branch 7) is still the canonical install path.

cran_pkgs <- c("yaml")
bioc_pkgs <- c("SummarizedExperiment", "recount3")

install_if_missing <- function(pkgs, repo = c("CRAN", "Bioc")) {
  repo <- match.arg(repo)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible())
  message("Installing ", repo, " packages: ", paste(missing, collapse = ", "))
  if (repo == "CRAN") {
    install.packages(missing, repos = "https://cloud.r-project.org")
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing, update = FALSE, ask = FALSE)
  }
}

install_if_missing(cran_pkgs, "CRAN")
install_if_missing(bioc_pkgs, "Bioc")

suppressPackageStartupMessages({
  library(yaml)
  library(SummarizedExperiment)
  library(recount3)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
cfg <- yaml::read_yaml(config_path)

`%||%` <- function(a, b) if (is.null(a)) b else a

# All paths in config.yaml are relative to the repo root.
# Resolve repo root from config file location and setwd so relative paths work
# regardless of whether the script is run from workflow/scripts/ or repo root.
config_abs <- normalizePath(config_path, mustWork = TRUE)
repo_root  <- dirname(dirname(config_abs))   # <repo>/config/config.yaml -> <repo>
setwd(repo_root)
message("Working directory set to repo root: ", repo_root)

set.seed(cfg$runtime$seed)
dir.create(dirname(cfg$paths$raw_counts), recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Pull or load counts ----------------------------------------------------

source_type <- cfg$dataset$source

if (source_type == "recount3") {
  projects <- recount3::available_projects()
  proj_info <- projects[
    projects$project == cfg$dataset$project &
      projects$file_source == "tcga" &
      projects$project_type == "data_sources",
  ]
  if (nrow(proj_info) != 1) {
    stop("recount3: expected 1 project row for ", cfg$dataset$project,
         " (tcga), got ", nrow(proj_info))
  }
  rse <- recount3::create_rse(proj_info)
  assay(rse, "counts") <- recount3::transform_counts(rse)

  counts   <- assay(rse, "counts")
  row_meta <- as.data.frame(rowData(rse))
  meta_df  <- as.data.frame(colData(rse))

  # Decide what to do with TCGA "Metastatic" samples before assigning condition.
  # Without an explicit policy they get silently lumped into the tumor arm.
  metastatic_action <- cfg$qc$metastatic %||% "drop"
  if (!metastatic_action %in% c("drop", "tumor", "separate")) {
    stop("qc.metastatic must be one of: drop, tumor, separate")
  }
  st        <- meta_df$tcga.cgc_sample_sample_type
  is_normal <- grepl("Normal",     st, ignore.case = TRUE)
  is_meta   <- grepl("Metastatic", st, ignore.case = TRUE)
  message("Sample types: ",
          paste(names(table(st)), table(st), sep = "=", collapse = ", "))

  if (metastatic_action == "drop" && any(is_meta)) {
    message("Metastatic handling: dropping ", sum(is_meta), " metastatic sample(s)")
    keep_s    <- !is_meta
    counts    <- counts[, keep_s, drop = FALSE]
    meta_df   <- meta_df[keep_s, , drop = FALSE]
    st        <- st[keep_s]
    is_normal <- is_normal[keep_s]
    is_meta   <- is_meta[keep_s]
  }

  condition <- ifelse(
    is_normal, cfg$dataset$reference_level,
    ifelse(is_meta & metastatic_action == "separate",
           "metastatic", cfg$dataset$treatment_level)
  )

  barcode <- meta_df$tcga.tcga_barcode
  samples <- data.frame(
    sample_id   = colnames(counts),
    barcode     = barcode,
    patient_id  = substr(barcode, 1, 12),       # TCGA-XX-YYYY
    condition   = condition,
    sample_type = st,
    stringsAsFactors = FALSE
  )

  # Pull configured clinical trait columns from recount3 colData. TCGA mixes
  # "[Not Available]" tokens into otherwise-numeric fields, so we coerce those
  # to NA and try numeric — downstream 03_wgcna.R then knows what's numeric vs
  # categorical without re-parsing strings.
  trait_map <- cfg$traits$recount3_columns
  if (length(trait_map) > 0) {
    na_tokens <- c("", "NA", "[Not Available]", "[Unknown]",
                   "[Not Evaluated]", "[Not Applicable]")
    clean_trait <- function(v) {
      if (is.factor(v)) v <- as.character(v)
      if (is.character(v)) {
        v[v %in% na_tokens] <- NA
        num <- suppressWarnings(as.numeric(v))
        n_parsed <- sum(!is.na(num))
        n_present <- sum(!is.na(v))
        if (n_present > 0 && n_parsed >= 0.5 * n_present) v <- num
      }
      v
    }
    for (short_name in names(trait_map)) {
      src_col <- trait_map[[short_name]]
      if (src_col %in% colnames(meta_df)) {
        samples[[short_name]] <- clean_trait(meta_df[[src_col]])
      } else {
        message("Trait '", short_name, "': source column '", src_col,
                "' not in recount3 colData; skipping")
      }
    }
  }
} else if (source_type == "local") {
  counts <- as.matrix(read.table(cfg$paths$raw_counts, header = TRUE,
                                 row.names = 1, sep = "\t", check.names = FALSE))
  samples <- read.table(cfg$paths$sample_metadata, header = TRUE,
                        sep = "\t", stringsAsFactors = FALSE)
  # Don't assume the metadata file is in counts-column order.
  samples <- samples[match(colnames(counts), samples$sample_id), , drop = FALSE]
  stopifnot(!anyNA(samples$sample_id),
            identical(samples$sample_id, colnames(counts)))
  row_meta <- NULL
} else {
  stop("Unsupported dataset.source: ", source_type)
}

stopifnot(identical(samples$sample_id, colnames(counts)))
message("Loaded ", ncol(counts), " samples, ", nrow(counts), " features")
message("Condition table: ",
        paste(names(table(samples$condition)),
              table(samples$condition), sep = "=", collapse = ", "))

# ---- De-duplicate replicate aliquots (TCGA only) ----------------------------
# recount3 sometimes carries multiple aliquots per (patient, sample_type). For
# DE this is usually harmless but for WGCNA it inflates apparent co-expression.
# Sort barcodes lexicographically and keep the first aliquot per (patient, condition).

if (isTRUE(cfg$qc$dedupe_aliquots) && !is.null(samples$patient_id)) {
  ord     <- order(samples$barcode)
  samples <- samples[ord, , drop = FALSE]
  counts  <- counts[, ord, drop = FALSE]
  dup     <- duplicated(paste(samples$patient_id, samples$condition, sep = "|"))
  if (any(dup)) {
    message("Aliquot dedup: dropped ", sum(dup), " duplicate aliquot(s)")
    samples <- samples[!dup, , drop = FALSE]
    counts  <- counts[, !dup, drop = FALSE]
  } else {
    message("Aliquot dedup: no duplicates found")
  }
}

# ---- Sample QC: library size flag -------------------------------------------

min_lib   <- cfg$qc$min_library_size %||% 0
libsize   <- colSums(counts)
samples$lib_size     <- libsize
samples$lib_size_low <- libsize < min_lib
message("Library size QC: ", sum(samples$lib_size_low),
        " sample(s) below threshold ", min_lib,
        " (min=", min(libsize), ", median=", round(median(libsize)),
        ", max=", max(libsize), ")")

# ---- Gene biotype filter ----------------------------------------------------
# Drops snoRNAs, lncRNAs, pseudogenes, etc. — only the protein-coding set is
# carried forward. Skipped if row_meta is unavailable (local counts table).

if (!is.null(row_meta) && isTRUE(cfg$prefilter$protein_coding_only)) {
  if (!"gene_type" %in% colnames(row_meta)) {
    warning("protein_coding_only=TRUE but rowData has no gene_type column; skipping")
  } else {
    pc <- row_meta$gene_type == "protein_coding"
    message("Biotype filter: kept ", sum(pc), " protein-coding / ",
            length(pc), " features")
    counts   <- counts[pc, , drop = FALSE]
    row_meta <- row_meta[pc, , drop = FALSE]
  }
}

# ---- Pre-filter: drop low-expression genes (CPM-based) ----------------------
# CPM normalises out library-size variation, so a gene with 10 counts in a
# shallow sample isn't treated the same as 10 counts in a deeply-sequenced one.

cpm  <- t(t(counts) / libsize) * 1e6
keep <- rowSums(cpm >= cfg$prefilter$min_cpm) >= cfg$prefilter$min_samples
message("Low-expression filter (CPM >= ", cfg$prefilter$min_cpm,
        " in >= ", cfg$prefilter$min_samples, " samples): kept ",
        sum(keep), " / ", length(keep), " features")
counts <- counts[keep, , drop = FALSE]
if (!is.null(row_meta)) row_meta <- row_meta[keep, , drop = FALSE]

stopifnot(identical(samples$sample_id, colnames(counts)))

# ---- Write deliverables -----------------------------------------------------

write.table(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
            file = cfg$paths$raw_counts,
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(samples,
            file = cfg$paths$sample_metadata,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", cfg$paths$raw_counts)
message("Wrote ", cfg$paths$sample_metadata)

if (!is.null(row_meta)) {
  genes_path <- file.path(dirname(cfg$paths$raw_counts), "genes.tsv")
  gene_cols  <- intersect(c("gene_name", "gene_type", "bp_length"), colnames(row_meta))
  genes_out  <- data.frame(gene_id = rownames(counts),
                           row_meta[, gene_cols, drop = FALSE],
                           stringsAsFactors = FALSE)
  write.table(genes_out, file = genes_path,
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", genes_path)
}

saveRDS(list(counts = counts, samples = samples, row_meta = row_meta),
        file = file.path(cfg$paths$rds_dir, "load_qc.rds"))
