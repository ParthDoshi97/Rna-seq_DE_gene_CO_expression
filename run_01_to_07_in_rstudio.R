# Run the analysis scripts from RStudio.
#
# The numbered scripts accept a config path when called with Rscript. When they
# are sourced interactively, this helper runs each script from its own folder so
# the default "../../config/config.yaml" path resolves correctly.

scripts <- file.path(
  "workflow",
  "scripts",
  sprintf("%02d_%s.R", 1:7, c(
    "load_and_qc",
    "differential_expression",
    "wgcna",
    "enrichment",
    "module_enrichment",
    "module_preservation",
    "hub_genes"
  ))
)

missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts) > 0) {
  stop("Missing script(s):\n  ", paste(missing_scripts, collapse = "\n  "))
}

for (script in scripts) {
  message("\n>>> Running ", script)
  source(script, chdir = TRUE)
}

message("\nAnalysis scripts 01-07 completed.")
