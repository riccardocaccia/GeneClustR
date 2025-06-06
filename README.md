# GeneClustR

**Project 2: R Package for Gene Expression Analysis and Clustering**

![R](https://img.shields.io/badge/R-≥4.1.0-blue?logo=R)
![Bioconductor](https://img.shields.io/badge/Bioconductor-compatible-green?logo=bioconductor)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

## Overview

**GeneClustR** is an R package designed for the analysis of gene expression data and the clustering of genes based on their expression profiles across different conditions or time points.

It supports multiple input formats and offers a full pipeline from preprocessing to clustering, visualization, and enrichment analysis — ideal for researchers exploring high-throughput transcriptomic datasets.

---

## Features

- **Data Import & Preprocessing**
  - Load expression data from `.csv` or Bioconductor `ExpressionSet` objects.
  - Normalize using log-transformation, quantile normalization, and more.

- **Exploratory Data Analysis**
  - Visualize data distributions.
  - Identify and optionally remove outlier genes.

- **Clustering Algorithms**
  - Hierarchical clustering
  - K-means
  - (Optionally extendable to other methods)

- **Visualization Tools**
  - Heatmaps of expression patterns
  - PCA plots
  - Dendrograms

- **Functional Enrichment Analysis**
  - Gene Set Enrichment Analysis (GSEA) on clusters
  - Integration with pathway databases

- **Documentation**
  - Comprehensive help files (`?function_name`)
  - Example datasets included
  - Vignettes for full pipeline demonstrations

---

