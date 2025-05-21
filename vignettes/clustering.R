## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(GeneClustR)
library(ggplot2)
library(plotly)
library(cluster)
library(factoextra)

## -----------------------------------------------------------------------------
set.seed(123)
mat <- matrix(rnorm(300), nrow = 30, ncol = 10)
rownames(mat) <- paste0("Gene", 1:30)
colnames(mat) <- paste0("Sample", 1:10)

# Preprocessing
mat_prep <- preprocess(mat, normalize = TRUE, method = "log")

# Rimuovi righe con valori non finiti (Inf, -Inf, NA, NaN)
mat_prep <- mat_prep[apply(mat_prep, 1, function(x) all(is.finite(x))), ]

## -----------------------------------------------------------------------------
cluster_number(mat_prep, method = "elbow", max_clusters = 10)


## -----------------------------------------------------------------------------
cluster_number(mat_prep, method = "gapstats", max_clusters = 10)


## -----------------------------------------------------------------------------
kmeans_res <- clustering(mat_prep, method = "kmeans", centers = 3)
kmeans_res$plot


## -----------------------------------------------------------------------------
kmedoids_res <- clustering(mat_prep, method = "kmedoids", centers = 3, metric = "euclidean")
kmedoids_res$plot


## -----------------------------------------------------------------------------
hclust_res <- clustering(mat_prep, method = "hierarchical", centers = 3)
hclust_res$plot  # PCA view
hclust_res$dendrogram  # Dendrogramma


## -----------------------------------------------------------------------------
# eps e minpts vanno trovati tramite sperimentazione (es. kNNdistplot)
dbscan_res <- clustering(mat_prep, method = "DBSCAN", eps = 1, minpts = 2)
dbscan_res$plot


