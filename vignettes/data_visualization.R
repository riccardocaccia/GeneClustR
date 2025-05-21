## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(GeneClustR)

## -----------------------------------------------------------------------------
# Simuliamo una matrice di espressione genica
set.seed(123)
data_matrix <- matrix(rnorm(200), nrow = 20, ncol = 10)
rownames(data_matrix) <- paste0("Gene", 1:20)


## -----------------------------------------------------------------------------
# Eseguiamo il clustering con metodo kmeans e 3 cluster
clustering_result <- clustering(data_matrix, method = "kmeans", centers = 3)

## -----------------------------------------------------------------------------
# 1. Visualizzazione del plot interattivo (PCA)
clustering_result$plot


## -----------------------------------------------------------------------------
# 2. Visualizzazione PCA (dati)
head(clustering_result$pca_result)


## -----------------------------------------------------------------------------
# 3. Riepilogo delle dimensioni dei cluster
clustering_result$cluster_information


## -----------------------------------------------------------------------------
# 4. Dendrogramma (solo per hierarchical clustering)
# clustering_result$dendrogram

## -----------------------------------------------------------------------------
# 5. Assegnazione dei cluster
clustering_result$clusters

