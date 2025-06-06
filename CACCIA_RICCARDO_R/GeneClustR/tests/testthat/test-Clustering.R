test_that("cluster_number returns elbow plot", {
  set.seed(1)
  mat <- matrix(rnorm(100), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("sample", 1:10)

  expect_s3_class(cluster_number(mat, method = "elbow"), "gg")
})

test_that("cluster_number returns gapstats plot", {
  skip_on_cran()  # Optional: skip long test on CRAN
  set.seed(1)
  mat <- matrix(rnorm(100), nrow = 10)
  expect_s3_class(cluster_number(mat, method = "gapstats"), "gg")
})

test_that("cluster_number throws error on invalid method", {
  mat <- matrix(rnorm(100), nrow = 10)
  expect_error(cluster_number(mat, method = "invalid"))
})


test_that("kmeans clustering returns expected output", {
  set.seed(42)
  mat <- matrix(rnorm(100), nrow = 10)
  res <- clustering(mat, method = "kmeans", centers = 3)
  expect_type(res, "list")
  expect_s3_class(res$plot, "plotly")
  expect_true("pca_result" %in% names(res))
})

test_that("kmedoids clustering works with metric", {
  set.seed(42)
  mat <- matrix(rnorm(100), nrow = 10)
  res <- clustering(mat, method = "kmedoids", centers = 3, metric = "euclidean")
  expect_type(res, "list")
  expect_s3_class(res$plot, "plotly")
})

test_that("hierarchical clustering returns dendrogram", {
  set.seed(42)
  mat <- matrix(rnorm(100), nrow = 10)
  res <- clustering(mat, method = "hierarchical", centers = 3)
  expect_s3_class(res$dendrogram, "gg")
  expect_s3_class(res$plot, "plotly")
})

test_that("DBSCAN clustering works with eps and minpts", {
  set.seed(42)
  mat <- matrix(rnorm(100), nrow = 10)
  res <- clustering(mat, method = "DBSCAN", eps = 1, minpts = 2)
  expect_type(res, "list")
  expect_s3_class(res$plot, "plotly")
})

test_that("clustering stops on invalid method", {
  mat <- matrix(rnorm(100), nrow = 10)
  expect_error(clustering(mat, method = "invalid"))
})
