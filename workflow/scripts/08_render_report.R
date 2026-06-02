#!/usr/bin/env Rscript
# 08_render_report.R
# Render the Quarto report to docs/index.html so GitHub Pages can serve it.

suppressPackageStartupMessages({
  library(yaml)
  library(quarto)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
if (!file.exists(config_path) && file.exists(file.path("config", "config.yaml"))) {
  config_path <- file.path("config", "config.yaml")
}
config_abs <- normalizePath(config_path, mustWork = TRUE)
cfg <- yaml::read_yaml(config_abs)

# All paths in config.yaml are relative to the repo root.
# Resolve repo root from config file location and setwd so relative paths work
# regardless of whether the script is run from workflow/scripts/ or repo root.
repo_root  <- dirname(dirname(config_abs))   # <repo>/config/config.yaml -> <repo>
setwd(repo_root)
message("Working directory set to repo root: ", repo_root)

required_outputs <- c(
  cfg$paths$de_results_csv,
  cfg$paths$module_assignments_csv,
  cfg$paths$enrichment_go_csv,
  cfg$paths$enrichment_kegg_csv,
  cfg$paths$gsea_go_csv,
  cfg$paths$gsea_kegg_csv,
  cfg$paths$module_annotation_summary_csv,
  cfg$paths$module_preservation_tsv,
  cfg$paths$hub_genes_ranked_csv,
  cfg$paths$hub_genes_priority_csv
)
required_outputs <- required_outputs[!vapply(required_outputs, is.null, logical(1))]
missing_outputs <- required_outputs[!file.exists(required_outputs)]
if (length(missing_outputs) > 0) {
  warning("Some expected report inputs are missing:\n  ",
          paste(missing_outputs, collapse = "\n  "))
}

dir.create(dirname(cfg$paths$report_html), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cfg$paths$report_qmd), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cfg$paths$report_qmd)) {
  stop("Report template not found: ", cfg$paths$report_qmd)
}

quarto::quarto_render(
  input = cfg$paths$report_qmd,
  output_format = "html",
  execute_params = list(config_path = config_abs)
)

# Quarto renders next to the .qmd by default; move to docs/index.html.
src_html <- sub("\\.qmd$", ".html", cfg$paths$report_qmd)
if (src_html != cfg$paths$report_html) {
  file.copy(src_html, cfg$paths$report_html, overwrite = TRUE)
  file.remove(src_html)
}
message("Rendered ", cfg$paths$report_html)
