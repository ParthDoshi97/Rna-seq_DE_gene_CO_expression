#!/usr/bin/env Rscript
# 05_render_report.R
# Render the Quarto report to docs/index.html so GitHub Pages can serve it.
# The .qmd template is added in branch 6 (feat/report).

suppressPackageStartupMessages({
  library(yaml)
  library(quarto)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "../../config/config.yaml"
cfg <- yaml::read_yaml(config_path)

dir.create(dirname(cfg$paths$report_html), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cfg$paths$report_qmd)) {
  stop("Report template not found: ", cfg$paths$report_qmd,
       " (added in branch feat/report)")
}

quarto::quarto_render(
  input = cfg$paths$report_qmd,
  output_format = "html",
  execute_params = list(config_path = normalizePath(config_path))
)

# Quarto renders next to the .qmd by default; move to docs/index.html.
src_html <- sub("\\.qmd$", ".html", cfg$paths$report_qmd)
if (src_html != cfg$paths$report_html) {
  file.copy(src_html, cfg$paths$report_html, overwrite = TRUE)
  file.remove(src_html)
}
message("Rendered ", cfg$paths$report_html)
