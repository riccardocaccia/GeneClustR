## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(GeneClustR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(ggplot2)

## -----------------------------------------------------------------------------
set.seed(123)

# Simuliamo dati di espressione per 100 geni e 10 campioni
genes <- paste0("ENSG", sprintf("%011d", 1:100))
samples <- paste0("Sample", 1:10)
mat <- matrix(rnorm(1000), nrow = 100, dimnames = list(genes, samples))

# Preprocessamento semplice
mat <- preprocess(mat, normalize = TRUE, method = "log")

# Etichette di gruppo
groups <- rep(c("A", "B"), each = 5)

## -----------------------------------------------------------------------------
gsea_res_no <- gsea_analysis(mat, keyType = "ENSEMBL", condition = "no")
gsea_res_no$plot


## -----------------------------------------------------------------------------
gsea_res_cond <- gsea_analysis(mat, keyType = "ENSEMBL", condition = "yes", groups = groups)

if (!is.null(gsea_res_cond$plot)) {
  gsea_res_cond$plot
} else {
  cat("Nessun termine GO significativamente arricchito trovato.")
}


## -----------------------------------------------------------------------------
# Visualizza le prime righe dei risultati
head(gsea_res_cond$gsea@result)


