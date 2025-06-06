# Test 1: Simple clustering (KMeans)
test_that("data_visualization returns correct output for KMeans", {
  # Eseguiamo un clustering di esempio
  mat <- matrix(rnorm(20 * 5), nrow = 5)  # 5 features, 20 samples (after transpose)
  clustering_result <- GeneClustR::clustering(mat, method = "kmeans", centers = 3)

  # Test per il plot
  expect_s3_class(data_visualization(clustering_result, what = "plot"), "plotly")  # Cambiato a plotly

  # Test PCA
  expect_s3_class(data_visualization(clustering_result, what = "pca"), "data.frame")

  #cluster summary
  expect_type(data_visualization(clustering_result, what = "summary"), "list")

  # Test clusters
  expect_type(data_visualization(clustering_result, what = "clusters"), "integer")
})

# Test 2: Clustering Hierarchical
test_that("data_visualization returns correct output for Hierarchical Clustering", {
  # Hierarchical clus example
  mat <- matrix(rnorm(20 * 5), nrow = 5)  # 5 features, 20 samples
  clustering_result <- GeneClustR::clustering(mat, method = "hierarchical", centers = 3)

  # dendrogram
  expect_s3_class(data_visualization(clustering_result, what = "dendrogram"), "gg")

  # Test PCA
  expect_s3_class(data_visualization(clustering_result, what = "pca"), "data.frame")

  # Test cluster summary
  expect_type(data_visualization(clustering_result, what = "summary"), "list")
})

# Test 3: DBSCAN
test_that("data_visualization returns correct output for DBSCAN", {
  # Eseguiamo un clustering DBSCAN di esempio
  mat <- matrix(rnorm(20 * 5), nrow = 5)  # 5 features, 20 samples
  clustering_result <- GeneClustR::clustering(mat, method = "DBSCAN", eps = 0.5, minpts = 5)

  # Test DBSCAN
  expect_type(data_visualization(clustering_result, what = "dbscan"), "list")
})
