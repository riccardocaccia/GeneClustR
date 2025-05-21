test_that("CorrHeatmap analysis returns summary and a plotly object", {
  mat <- matrix(rnorm(100), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("sample", 1:10)

  res <- explorative_analysis(mat, analysis_type = "CorrHeatmap")
  expect_type(res, "list")
  expect_true("summary" %in% names(res))
  expect_s3_class(res$plot, "plotly")
})

test_that("scatter analysis requires valid genes", {
  mat <- matrix(rnorm(100), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("sample", 1:10)

  expect_error(explorative_analysis(mat, analysis_type = "scatter"))
  expect_error(explorative_analysis(mat, analysis_type = "scatter", gene1 = "gene1", gene2 = "NOT_IN_DATA"))
  expect_s3_class(explorative_analysis(mat, analysis_type = "scatter", gene1 = "gene1", gene2 = "gene2"), "ggplot")
})

test_that("lm analysis stops with invalid genes", {
  mat <- matrix(rnorm(100), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("sample", 1:10)

  expect_error(explorative_analysis(mat, analysis_type = "lm", gene1 = "gene1"))  # solo un gene
  expect_error(explorative_analysis(mat, analysis_type = "lm", gene1 = "gene1", gene2 = "nope"))
})

test_that("PCA returns expected structure", {
  mat <- matrix(rnorm(100), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("sample", 1:10)

  res <- explorative_analysis(mat, analysis_type = "pca")
  expect_type(res, "list")
  expect_s3_class(res$plot, "ggplot")
  expect_true("summary" %in% names(res))
  expect_true("pca_df" %in% names(res))
  expect_true("relevant_genes" %in% names(res))
})

test_that("invalid analysis_type triggers error", {
  mat <- matrix(rnorm(100), nrow = 10)
  expect_error(explorative_analysis(mat, analysis_type = "invalid"))
})
