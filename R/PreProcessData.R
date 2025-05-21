#' Preprocess DNA expression data
#'
#' Given genes expression data (as matrix or ExpressionSet), preprocess giving a normalized df
#'
#' @param data Matrix or ExpressionSet, of expression values
#' @param normalize Logical, decide if you want to normalize your data, default=TRUE
#' @param method Character, choose the method to normalize between ('log', 'quantile'), default='log'
#' @param batch_effect Logical, decide whether to apply batch effect correction, default=FALSE
#' @param batch Vector, indicating the batch assignment, not required for ExpressionSet
#'
#' @importFrom limma normalizeQuantiles
#' @importFrom sva ComBat
#' @importFrom Biobase exprs pData
#' @return A normalized matrix or ExpressionSet
#' @export
preprocess <- function(data, normalize = TRUE, method = 'log', batch_effect = FALSE, batch = NULL) {

  # ExpressionSet handling
  if (inherits(data, 'ExpressionSet')) {
    expression <- Biobase::exprs(data)

    if (normalize) {
      if (method == 'log') {
        expression <- log2(expression + 1)
      } else if (method == 'quantile') {
        expression <- normalizeQuantiles(expression)
      } else {
        stop('Provide "log" or "quantile" as normalization method.')
      }
    }

    if (batch_effect) {
      if (!"batch" %in% colnames(Biobase::pData(data))) {
        stop("Batch effect correction requested, but no 'batch' column in pData.")
      }
      batch <- Biobase::pData(data)$batch
      expression <- sva::ComBat(dat = expression, batch = batch)
    }

    Biobase::exprs(data) <- expression
    return(data)
  }

  # Matrix handling
  else if (is.matrix(data)) {

    # Controlla se ci sono valori non numerici (NA, NaN)
    if (!all(is.finite(data))) {
      warning("I dati contengono valori non finiti (NA, NaN). Saranno sostituiti con la media della riga.")
    }

    # Replace NA/NaN with row means
    for (i in 1:nrow(data)) {
      na_value <- is.na(data[i, ]) | is.nan(data[i, ])
      if (any(na_value)) {
        row_mean <- mean(data[i, !na_value], na.rm = TRUE)
        data[i, na_value] <- row_mean
      }
    }

    if (normalize) {
      if (method == 'log') {
        data[data <= 0] <- 1e-3  # Assicura che tutti i valori siano > 0
        data <- log2(data)
      } else if (method == 'quantile') {
        data <- normalizeQuantiles(data)
      } else {
        stop('Provide "log" or "quantile" as normalization method.')
      }

      # Rimuovi righe con valori non finiti o varianza zero
      data <- data[apply(data, 1, function(x) all(is.finite(x)) && sd(x) > 0), ]

      if (nrow(data) == 0 || ncol(data) == 0) {
        stop("La matrice risultante è vuota dopo il preprocessing. Controlla i dati iniziali.")
      }
    }

    if (batch_effect) {
      if (is.null(batch)) stop("Provide a batch vector to apply batch effect correction.")
      data <- sva::ComBat(dat = data, batch = batch)
    }

    return(data)
  }

  else {
    stop("Unsupported data type. Must be a matrix or ExpressionSet.")
  }
}
