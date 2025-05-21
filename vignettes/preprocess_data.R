## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(GeneClustR)
library(Biobase)
library(sva)

## -----------------------------------------------------------------------------
# Simuliamo una matrice di espressione
set.seed(42)
data_matrix <- matrix(rnorm(1000, mean = 100, sd = 20), nrow = 100, ncol = 10)
data_matrix[data_matrix <= 0] <- 1e-3
rownames(data_matrix) <- paste0("Gene", 1:100)

# Preprocessing log
preprocessed_log <- preprocess(data_matrix, normalize = TRUE, method = "log")

# Preprocessing quantile
preprocessed_quantile <- preprocess(data_matrix, normalize = TRUE, method = "quantile")

# Visualizzazione
head(preprocessed_log[, 1:5])


## -----------------------------------------------------------------------------
# Simuliamo una matrice di espressione
set.seed(42)
data_matrix <- matrix(rnorm(1000, mean = 100, sd = 20), nrow = 100, ncol = 10)
data_matrix[data_matrix <= 0] <- 1e-3
rownames(data_matrix) <- paste0("Gene", 1:100)

# Preprocessing log
preprocessed_log <- preprocess(data_matrix, normalize = TRUE, method = "log")

# Preprocessing quantile
preprocessed_quantile <- preprocess(data_matrix, normalize = TRUE, method = "quantile")

# Visualizzazione
head(preprocessed_log[, 1:5])


## -----------------------------------------------------------------------------
# Simuliamo dati con batch
batch_vector <- rep(c("A", "B"), length.out = ncol(data_matrix))

# Preprocessing con batch effect
preprocessed_batch <- preprocess(data_matrix, normalize = TRUE, method = "log",
                                 batch_effect = TRUE, batch = batch_vector)


