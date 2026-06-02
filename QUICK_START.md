# Quick Start

This guide shows how to rerun the TCGA-BRCA RNA-seq analysis from an RStudio
session or from the command line.

## 1. Start From The Repo Root

```r
setwd("path/to/Rna-seq_DE_gene_CO_expression")
```

The main configuration file is:

```text
config/config.yaml
```

All scripts resolve paths relative to the repo root.

## 2. RStudio Run

For the analysis steps that generate the main results:

```r
source("run_01_to_07_in_rstudio.R")
```

Then render the report:

```r
source("workflow/scripts/08_render_report.R")
```

The final HTML report is written to:

```text
docs/index.html
```

When browsing the repository on GitHub, opening `docs/index.html` shows the
HTML source. To view the rendered webpage, merge the report into `main`, enable
GitHub Pages from the `main` branch and `/docs` folder, then use:

```text
https://parthdoshi97.github.io/Rna-seq_DE_gene_CO_expression/
```

## 3. Script Order

| Step | Script | Main output |
|---|---|---|
| 01 | `01_load_and_qc.R` | `data/counts.tsv`, `data/samples.tsv`, `data/genes.tsv` |
| 02 | `02_differential_expression.R` | `results/tables/de_results.tsv` |
| 03 | `03_wgcna.R` | `results/tables/module_assignments.tsv` |
| 04 | `04_enrichment.R` | DE GO/KEGG ORA and GSEA tables/figures |
| 05 | `05_module_enrichment.R` | Per-module enrichment tables/figures |
| 06 | `06_module_preservation.R` | `module_preservation.tsv` |
| 07 | `07_hub_genes.R` | ranked and priority hub-gene tables |
| 08 | `08_render_report.R` | `docs/index.html` |

## 4. Command-line Run

If `Rscript` is available on your PATH:

```bash
Rscript workflow/scripts/01_load_and_qc.R config/config.yaml
Rscript workflow/scripts/02_differential_expression.R config/config.yaml
Rscript workflow/scripts/03_wgcna.R config/config.yaml
Rscript workflow/scripts/04_enrichment.R config/config.yaml
Rscript workflow/scripts/05_module_enrichment.R config/config.yaml
Rscript workflow/scripts/06_module_preservation.R config/config.yaml
Rscript workflow/scripts/07_hub_genes.R config/config.yaml
Rscript workflow/scripts/08_render_report.R config/config.yaml
```

## 5. Expected Key Results

After a successful run, check these files:

```text
results/tables/de_results.tsv
results/tables/module_annotation_summary.tsv
results/tables/module_preservation.tsv
results/tables/hub_genes_priority_modules.tsv
docs/index.html
```

Important current interpretation:

- greenyellow is **preserved** (`Zsummary = 8.52`), not a not-preserved module.
- greenyellow is prioritized because of strong functional enrichment in
  synaptic/neurotransmitter signaling.
- The hub table records this as `priority_reason = enriched`.

## 6. Common Issues

### `Rscript` is not found

Run the scripts from RStudio, or call the full Rscript path, for example:

```powershell
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" workflow/scripts/08_render_report.R config/config.yaml
```

### Report renders but plots are missing

Rerun the current `08_render_report.R`. The report template resolves figure
paths from the repo root, so plots should be embedded in `docs/index.html`.

### Module preservation is slow

The preservation step uses WGCNA permutations. For a quick test, temporarily
lower this value in `config/config.yaml`:

```yaml
preservation:
  n_permutations: 50
```

For final analysis results, keep the configured value higher.
