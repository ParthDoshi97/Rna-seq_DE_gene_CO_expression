# 01_load_and_qc.R
# Pull TCGA-BRCA counts via recount3 (or load a local counts table),
# build the sample table, apply low-count pre-filtering, and emit
# a counts.tsv + samples.tsv pair that downstream scripts consume.

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

  counts <- assay(rse, "counts")
  meta_df <- as.data.frame(colData(rse))

  # TCGA sample_type comes from cgc_sample_sample_type
  meta_df$condition <- ifelse(
    grepl("Normal", meta_df$tcga.cgc_sample_sample_type, ignore.case = TRUE),
    cfg$dataset$reference_level,
    cfg$dataset$treatment_level
  )

  samples <- data.frame(
    sample_id  = colnames(counts),
    condition  = meta_df$condition,
    sample_type = meta_df$tcga.cgc_sample_sample_type,
    stringsAsFactors = FALSE
  )
} else if (source_type == "local") {
  counts <- as.matrix(read.table(cfg$paths$raw_counts, header = TRUE,
                                 row.names = 1, sep = "\t", check.names = FALSE))
  samples <- read.table(cfg$paths$sample_metadata, header = TRUE,
                        sep = "\t", stringsAsFactors = FALSE)
} else {
  stop("Unsupported dataset.source: ", source_type)
}

stopifnot(all(samples$sample_id == colnames(counts)))
message("Loaded ", ncol(counts), " samples, ", nrow(counts), " genes")
message("Condition table: ",
        paste(names(table(samples$condition)),
              table(samples$condition), sep = "=", collapse = ", "))

# ---- Pre-filter: drop genes with too few counts -----------------------------

keep <- rowSums(counts >= cfg$prefilter$min_count) >= cfg$prefilter$min_samples
message("Pre-filter: kept ", sum(keep), " / ", length(keep), " genes")
counts <- counts[keep, , drop = FALSE]

# ---- Write deliverables -----------------------------------------------------

write.table(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
            file = cfg$paths$raw_counts,
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(samples,
            file = cfg$paths$sample_metadata,
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(list(counts = counts, samples = samples),
        file = file.path(cfg$paths$rds_dir, "load_qc.rds"))

message("Wrote ", cfg$paths$raw_counts)
message("Wrote ", cfg$paths$sample_metadata)
