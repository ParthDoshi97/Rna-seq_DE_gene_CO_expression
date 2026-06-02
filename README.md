# TCGA-BRCA RNA-seq Differential Expression and WGCNA Analysis Project

This repository is an end-to-end bulk RNA-seq analysis of TCGA breast invasive
carcinoma (BRCA) tumor and normal samples. The project combines DESeq2,
functional enrichment, WGCNA co-expression modules, module preservation, and
hub-gene prioritization into a reproducible analysis report.

The biological question is:

> Which tumor-associated expression programs and co-expression modules are
> visible in TCGA-BRCA, and which hub genes best represent those programs?

The repository is designed as a reproducible bioinformatics analysis that can
be inspected, rerun, and adapted. It emphasizes not only code that runs, but
also data provenance, interpretation, reproducibility, and limitations.

## Quick Links

- GitHub-readable report: [REPORT.md](REPORT.md)
- Report HTML source: [docs/index.html](docs/index.html)
- Live report: available after this branch is merged and GitHub Pages is enabled
- Pipeline config: [config/config.yaml](config/config.yaml)
- Run helper for RStudio: [run_01_to_07_in_rstudio.R](run_01_to_07_in_rstudio.R)
- Report renderer: [workflow/scripts/08_render_report.R](workflow/scripts/08_render_report.R)
- Data provenance notes: [data/README.md](data/README.md)
- Results guide: [results/README.md](results/README.md)

## Project Snapshot

| Analysis area | Current result |
|---|---:|
| TCGA-BRCA samples | 1,207 |
| Tumor samples | 1,095 |
| Normal samples | 112 |
| Genes tested for differential expression | 18,158 |
| Significant DE genes, padj < 0.05 and abs(log2FC) > 1 | 4,823 |
| Up-regulated genes | 2,854 |
| Down-regulated genes | 1,969 |
| Genes used for WGCNA | 5,000 |
| Non-grey WGCNA modules | 12 |
| Priority hub-gene rows | 188 |

## Main Biological Findings

The analysis recovers several expected and biologically interpretable breast
cancer programs:

- **Turquoise module:** chromosome segregation and cell-cycle biology; strongly
  associated with ER status. Top hub genes include `FOXA1`, `MLPH`, `ESR1`,
  `XBP1`, and `GATA3`.
- **Brown module:** lymphocyte activation and adaptive immune response.
- **Blue and pink modules:** extracellular matrix and collagen-related biology.
- **Black module:** vasculature development and angiogenesis.
- **Greenyellow module:** synaptic and neurotransmitter signaling,
  including hub genes such as `RUNDC3A`, `CPLX2`, `SYP`, `UNC13A`, `SNAP25`,
  and `SYT4`.

The greenyellow module is an important example of why the hub-prioritization
logic uses more than one gate. It is **preserved** in the normal comparison
network (`Zsummary = 8.52`), so it should not be described as tumor-specific
rewiring. However, it has strong coherent functional annotation
(`best_padj = 3.26e-11` for synaptic signaling terms), so it is retained as a
priority module through the enrichment gate.

## Analysis Workflow

The scripts are numbered in execution order.

| Step | Script | Purpose |
|---|---|---|
| 01 | `workflow/scripts/01_load_and_qc.R` | Download/load TCGA-BRCA data, filter genes, derive metadata |
| 02 | `workflow/scripts/02_differential_expression.R` | DESeq2 tumor-vs-normal differential expression |
| 03 | `workflow/scripts/03_wgcna.R` | Tumor-only WGCNA module detection and trait correlation |
| 04 | `workflow/scripts/04_enrichment.R` | GO/KEGG ORA and GSEA for DE results |
| 05 | `workflow/scripts/05_module_enrichment.R` | GO/KEGG enrichment for each WGCNA module |
| 06 | `workflow/scripts/06_module_preservation.R` | Tumor-vs-normal module preservation analysis |
| 07 | `workflow/scripts/07_hub_genes.R` | Hub-gene ranking and priority module selection |
| 08 | `workflow/scripts/08_render_report.R` | Render the Quarto HTML report |

## Reproducibility

The pipeline is controlled from a single configuration file:

```text
config/config.yaml
```

For an RStudio-oriented run, use:

```r
source("run_01_to_07_in_rstudio.R")
source("workflow/scripts/08_render_report.R")
```

The helper script installs missing R/Bioconductor packages into a project-local
`.R-lib` directory. A Conda-style dependency file is also provided:

```bash
conda env create -f environment.yml
conda activate rnaseq-deseq2-wgcna
```

Some steps are slow, especially module preservation, because they use WGCNA
permutations. The number of preservation permutations is configured in
`config/config.yaml`.

## Repository Structure

```text
config/                 # Pipeline paths, thresholds, and method settings
data/                   # Generated data files; see data/README.md
docs/                   # Rendered HTML report for GitHub Pages
report/                 # Quarto report source
results/                # Generated tables and figures; see results/README.md
workflow/scripts/       # Numbered R scripts for the analysis
run_01_to_07_in_rstudio.R
environment.yml
```

Large raw/intermediate data and generated result files are intentionally ignored
by git. The rendered report in `docs/index.html` is the compact shareable
analysis summary.

## Viewing The Reports

For immediate viewing in GitHub, open [REPORT.md](REPORT.md). This Markdown
version includes the main result tables and selected figures.

GitHub shows `docs/index.html` as source code when it is opened from the
repository file browser. To view it as a webpage, the report must be served by
GitHub Pages.

Before this branch is merged into `main`, the Pages URL may return `404`
because `docs/index.html` is only present on the PR branch. After the branch is
merged and Pages is enabled, the report URL will be:

```text
https://parthdoshi97.github.io/Rna-seq_DE_gene_CO_expression/
```

For a newly merged repository, enable Pages in GitHub with:

```text
Settings -> Pages -> Build and deployment -> Source: GitHub Actions
```

After Pages is enabled, `.github/workflows/deploy-pages.yml` deploys the
`docs/` folder, so `docs/index.html` is served as the site homepage.

## Method Choices

- **DESeq2** was used because RNA-seq counts are discrete and overdispersed.
  The workflow uses adjusted p-values and log2 fold-change thresholds to define
  significant DE genes.
- **WGCNA** was run on tumor samples so modules capture tumor expression
  structure rather than simply separating tumor from normal.
- **GO/KEGG enrichment** was used to translate gene lists and modules into
  biological themes.
- **Module preservation** compares tumor modules against normal samples to ask
  whether module structure is retained outside the tumor context.
- **Hub prioritization** keeps modules that are trait-associated,
  not-preserved, or strongly enriched. This prevents biologically coherent but
  modestly trait-correlated modules from disappearing.

## Limitations

- This is a public-data reanalysis, not a clinical validation study.
- TCGA bulk RNA-seq mixes tumor cells, stromal cells, immune cells, and other
  cell types. Module biology may reflect tumor microenvironment composition as
  well as cancer-cell-intrinsic programs.
- Clinical trait fields in TCGA can be sparse, noisy, or inconsistently coded.
- The analysis does not currently model all possible confounders such as batch,
  purity, treatment history, or sequencing center.
- Hub genes are candidates for interpretation, not proven drivers.

## Project Strengths

This project is structured to make the analysis easy to inspect and rerun:

- clear biological question and public data source
- numbered, reproducible scripts
- configuration-driven thresholds and paths
- QC, DE, enrichment, WGCNA, preservation, and hub-gene outputs
- rendered scientific report
- explicit interpretation and limitations

## License and Citation

This project is released under the Apache-2.0 license. Citation metadata is
provided in [CITATION.cff](CITATION.cff).
