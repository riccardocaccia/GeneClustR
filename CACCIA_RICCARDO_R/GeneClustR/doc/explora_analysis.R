## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(GeneClustR)
library(ggplot2)
library(plotly)

## -----------------------------------------------------------------------------
# Creazione matrice simulata
set.seed(10)
mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
rownames(mat) <- paste0("Gene", 1:20)
colnames(mat) <- paste0("Sample", 1:10)

# Preprocessing
mat_prep <- preprocess(mat, normalize = TRUE, method = "log")

# Cleaning: rimuove righe con NA/Inf e verifica che rimanga qualcosa
finite_rows <- apply(mat_prep, 1, function(x) all(is.finite(x)))
finite_rows[is.na(finite_rows)] <- FALSE  # evita NA

mat_prep <- mat_prep[finite_rows, , drop = FALSE]

# ⚠️ Check: la matrice non deve essere vuota
if (nrow(mat_prep) == 0 || ncol(mat_prep) == 0) {
  stop("La matrice pre-processata è vuota dopo la rimozione dei dati non finiti.")
}

mat_prep <- as.matrix(mat_prep)
options(repr.plot.width = 10, repr.plot.height = 8)

## -----------------------------------------------------------------------------
results <- explorative_analysis(mat_prep, analysis_type = "CorrHeatmap")
results$plot


## -----------------------------------------------------------------------------
pca_out <- explorative_analysis(mat_prep, analysis_type = "pca")
pca_out$plot

## -----------------------------------------------------------------------------
explorative_analysis(mat_prep, analysis_type = "scatter", gene1 = "Gene1", gene2 = "Gene2")


## -----------------------------------------------------------------------------
# Nota: questo produce output base R, non ggplot
explorative_analysis(mat_prep, analysis_type = "lm", gene1 = "Gene1", gene2 = "Gene2")


