# Benchmarking Confounder-aware Differential Abundance (xConDA) Methods in Microbiome Data

## Workflow
<!-- <img src="imgs/Study_design.png" width="600" alt="Study overview"> -->

## Introduction
We present a large-scale benchmark of **72 confounder-aware differential abundance analysis (DAA) strategies**, assembled from six input schemes (five normalization approaches plus raw counts) and 14 statistical models. 

## Software versions

The benchmark was run with the following software versions.

| Software | Version |
| --- | --- |
| R | 4.3.1 |
| Python | 3.8.18 |

## R package versions

The following versions correspond to R packages that are explicitly loaded or referenced in the benchmark scripts and were available in the recorded analysis environment.

| Package | Version |
| --- | --- |
| ALDEx2 | 1.34.0 |
| ANCOMBC | 2.4.0 |
| coin | 1.4-3 |
| compositions | 2.0-8 |
| corncob | 0.4.1 |
| DESeq2 | 1.42.1 |
| dplyr | 1.1.4 |
| edgeR | 4.0.16 |
| fastANCOM | 0.0.4 |
| ggplot2 | 3.5.2 |
| GUniFrac | 1.8 |
| limma | 3.58.1 |
| lmerTest | 3.1-3 |
| logging | 0.10.108 |
| Maaslin2 | 1.16.0 |
| metafor | 4.8.0 |
| metagenomeSeq | 1.43.0 |
| metap | 1.11 |
| MuMIn | 1.48.4 |
| paletteer | 1.6.0 |
| patchwork | 1.3.0 |
| pheatmap | 1.0.12 |
| pracma | 2.4.4 |
| purrr | 1.0.2 |
| SparseDOSSA2 | 0.99.2 |
| stringr | 1.5.1 |
| tibble | 3.2.1 |
| tidyverse | 2.0.0 |
| TreeSummarizedExperiment | 2.10.0 |
| VTwins | 0.1.0 |
<!-- 
- We systematically evaluate precision and sensitivity to identify **top performers** across 250 simulated scenarios.We then validate these strategies on three real-world datasets with approximate ground truth and assess cross-dataset consistency in 36 metagenomic datasets. We further probe robustness to sample size, prevalence, feature-effect magnitude, and the number, type, and strength of confounders.
- Beyond single-cohort evaluation, we further integrate the selected strategies with meta- or mega-analysis frameworks, evaluating both simulated batch-affected settings and real cohorts to recommend effective, **batch-robust DAA pipelines**.

We also provided **a DAA strategy benchmarking pipeline** that helps researchers choose an optimal strategy and integration framework for their datasets, and then execute the selected methods end to end to obtain robust differential microbes. 

* We provide a step-by-step **tutorial** with sample data for quick, reproducible use. 
* The pipeline is also available on the **xConDA** webserver.
-->
