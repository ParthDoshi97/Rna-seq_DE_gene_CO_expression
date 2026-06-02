# Data Provenance

This project analyzes public TCGA-BRCA RNA-seq data accessed through
`recount3`.

## Source

- Dataset source: TCGA via recount3
- Cancer type: Breast invasive carcinoma (BRCA)
- Genome/build context: GRCh38 / GENCODE v26 recount3 uniform reprocessing
- Contrast: tumor vs normal

The data-loading and QC script is:

```text
workflow/scripts/01_load_and_qc.R
```

The current analysis generated:

| File | Meaning |
|---|---|
| `counts.tsv` | Filtered count matrix used by downstream analysis |
| `samples.tsv` | Sample metadata, condition labels, and selected clinical traits |
| `genes.tsv` | Gene annotation used for gene symbols and biotypes |

These files are generated locally and are intentionally ignored by git because
they are large. Re-run script 01 to regenerate them.

## Data Handling Notes

- Only public TCGA/recount3 data are used.
- No private patient data were added to this repository.
- TCGA sample metadata can be incomplete or inconsistently coded, so clinical
  trait associations should be interpreted carefully.

