test_that("preprocess works correctly on a matrix", {
  # Crea una matrice finta di espressione
  set.seed(123)
  mat <- matrix(abs(rnorm(100)), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("sample", 1:10)

  # Test log normalization
  res_log <- preprocess(mat, normalize = TRUE, method = "log")
  expect_true(all(res_log >= 0))  # log2(x + 1) è sempre >= 0

  # Test quantile normalization
  res_quantile <- preprocess(mat, normalize = TRUE, method = "quantile")
  expect_equal(dim(res_quantile), dim(mat))
})

test_that("preprocess handles NA values by imputing row means", {
  mat <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2, byrow = TRUE)
  mat_res <- preprocess(mat, normalize = FALSE)
  expect_false(any(is.na(mat_res)))
})

test_that("preprocess errors on wrong normalization method", {
  mat <- matrix(abs(rnorm(100)), nrow = 10)
  expect_error(preprocess(mat, normalize = TRUE, method = "invalid_method"))
})

test_that("preprocess errors on missing batch vector", {
  mat <- matrix(abs(rnorm(100)), nrow = 10)
  expect_error(preprocess(mat, batch_effect = TRUE))
})
