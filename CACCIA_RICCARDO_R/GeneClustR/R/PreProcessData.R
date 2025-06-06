#' Preprocess DNA expression data
#'
#' Given genes expression data (as matrix), preprocess with different
#' normalization methods giving as output normalized data
#'
#' @param data Matrix, expression values
#' @param normalize Logical, decide if you want to normalize your data, default=TRUE
#' @param method Character, choose the method to normalize between ('log', 'quantile'), default='log'
#' @param batch_effect Logical, decide whether to apply batch effect correction, default=FALSE
#' @param batch Vector, indicating the batch assignment, required for batch effect correction
#' @param fallback_n Integer, number of random genes to retain if all get removed, default=5
#'
#' @importFrom limma normalizeQuantiles
#' @importFrom sva ComBat
#' @return A normalized matrix
#' @examples
#' # Simulate a matrix
#' set.seed(123)
#' expr_data <- matrix(abs(rnorm(600, mean = 10, sd = 3)), nrow = 100, ncol = 6)
#' colnames(expr_data) <- paste0("Sample", 1:6)
#'
#' # Default preprocessing with log2 normalization
#' norm_data <- preprocess(expr_data)
#'
#' # Preprocessing with quantile normalization
#' norm_data_q <- preprocess(expr_data, method = "quantile")
#'
#' # Preprocessing with batch effect (2 groups of 3 samples each)
#' batch <- c(1,1,1,2,2,2)
#' norm_batch_data <- preprocess(expr_data, batch_effect = TRUE, batch = batch)
#' @export
preprocess <- function(data, normalize = TRUE, method = 'log',
                       batch_effect = FALSE, batch = NULL, fallback_n = 5) {

  # Handle empty matrix
  handle_empty <- function(mat, fallback_n) {
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      message("The resulting matrix is empty,",
              fallback_n, " random rows will be selected from the original data")
      fallback_rows <- sample(1:nrow(data), fallback_n)
      return(data[fallback_rows, , drop = FALSE])
    }
    return(mat)
  }

  # Verify that data are in matrix format
  if (!is.matrix(data)) {
    stop("Input data must be in matrix format.")
  }

  # Hanfle NA o NaN values
  for (i in 1:nrow(data)) {
    na_value <- is.na(data[i, ]) | is.nan(data[i, ])
    if (any(na_value)) {
      row_mean <- mean(data[i, !na_value], na.rm = TRUE)
      data[i, na_value] <- row_mean
    }
  }

  # Normalization
  if (normalize) {
    if (method == 'log') {
      # Ensure that every value is >0 to apply log
      data[data <= 0] <- 1  # Change fallback value to 1 for log apllication
      data <- log2(data + 1)  # log2(x + 1) standard formula
    } else if (method == 'quantile') {
      data <- normalizeQuantiles(data)
    } else {
      stop('Fornisci "log" o "quantile" come metodo di normalizzazione.')
    }

    # Remove rows with null variance and Non finite values
    data <- data[apply(data, 1, function(x) all(is.finite(x)) && sd(x) > 0),
                 , drop = FALSE]

    # Handle empty matrix
    data <- handle_empty(data, fallback_n)
  }

  # Batch effect correction
  if (batch_effect) {
    if (is.null(batch)) {
      stop("Batch effect selected but 'batch vector' not specified.")
    }
    data <- sva::ComBat(dat = data, batch = batch)
  }

  return(data)
}
