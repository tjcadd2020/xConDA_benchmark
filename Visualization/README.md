# Main-figure visualization code

This directory contains panel-level R code corresponding approximately to the main figures in `Manuscript.docx`.

## Included panels

| Figure | Panels | Scripts |
| --- | --- | --- |
| Fig. 2 | a-e | Confounder adjustment, disease correlation, PERMANOVA variance partitioning, strategy ranking, and FDR-sensitivity trade-off |
| Fig. 3 | b-e | Oral reference recovery, oral and vaginal recovery scatter plots, and Crohn's disease stable recovery |
| Fig. 4 | a-d | Within-dataset consistency, dual-criterion prioritization, perturbation robustness, and challenging-setting performance |
| Fig. 5 | a and c | Integrated single-cohort performance profile and Crohn's disease case study |
| Fig. 6 | b, c, and e | Multi-cohort accuracy-consistency landscape, top-strategy profiles, and MMUPHin comparison |

## Data inputs

The scripts use benchmark result tables or saved R objects produced by the analysis workflow. Input filenames and required columns are documented at the beginning of each script. Place the corresponding inputs next to the script or pass an input path where the script supports command-line arguments.

Run a panel script from the repository root, for example:

```bash
Rscript Visualization/Fig2/Fig2b_disease_correlation.R
```

The repository's Docker image can be used to provide the R runtime and package environment described in the main README.
