#' Clustering number
#'
#' This function supports two techniques: elbow method and gap statistics. To
#' help estimate the optimal number of clusters to use in downstream clustering.
#'
#' @param data Matrix, cleaned and preprocessed gene expression or feature data
#' @param method Character, choose c(elbow, gapstats) estimation method
#' @param max_clusters Integer, maximum number of clusters to consider (default = 50)
#' @import factoextra
#' @import cluster
#'
#' @return a plot to understand how maany k you should use
#' @examples
#' # Simulated example
#' set.seed(42)
#' mat <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#'
#' # Elbow method
#' cluster_number(mat, method = "elbow")
#'
#' # Gap statistic
#' cluster_number(mat, method = "gapstats")
#' @export
cluster_number <- function(data, method = 'elbow', max_clusters = 50) {

  if (is.matrix(data)) {
    data <- as.data.frame(t(data))  # Transpose to have samples in rows
  } else {
    data <- as.data.frame(data)
  }

  n_samples <- nrow(data)
  K.max <- min(max_clusters, n_samples - 1)  # Ensure k < n

  if (K.max < 2) {
    stop("Not enough samples to compute clustering. Need at least 3 samples.")
  }

  # Elbow method
  if (method == 'elbow') {
    p <- fviz_nbclust(data, pam, method = "wss", k.max = K.max)
    return(p)
  }
  # GAp statsitsics method
  else if (method == 'gapstats') {
    #calculate gap statistic based on number of clusters
    gap_stat <- clusGap(data,
                        FUN = pam,
                        K.max = K.max, #max clusters to consider
                        B = 50) #total bootstrapped iterations

    #plot number of clusters vs. gap statistic
    p <- fviz_gap_stat(gap_stat)
    return(p)
  }

  # Invalid method
  else {
    stop('Invalid method')
  }
}
#' Clustering method
#'
#' This function performs clustering using different techniques
#' (k-means, k-medoids,hierarchical, DBSCAN) and returns outputs
#' including an interactive PCA plot,clustering assignments, and summary.
#'
#' @param data Matrix or data.frame, cleaned and normalized
#' @param method Character, method = c("kmeans", "hierarchical", "kmedoids", "DBSCAN")
#' @param centers Numerical, number of clusters
#' @param metric Character, distance metric for kmedoids ('euclidean','manhattan')
#' @param eps Numerical, epsilon value for DBSCAN
#' @param minpts Numerical, minimum points for DBSCAN
#' @import dbscan
#' @return A list with interactive plot, PCA result and clustering info
#' @examples
#' # Example using the iris dataset (excluding species label)
#' data(iris)
#' iris_data <- iris[, -5]  # Remove species label
#'
#' # K-means clustering
#' result_kmeans <- clustering(data = iris_data, method = "kmeans", centers = 3)
#' print(result_kmeans$cluster_information)
#'
#' # K-medoids clustering
#' result_kmedoids <- clustering(data = iris_data, method = "kmedoids", centers = 3, metric = "euclidean")
#' print(result_kmedoids$cluster_information)
#'
#' # Hierarchical clustering
#' result_hierarchical <- clustering(data = iris_data, method = "hierarchical", centers = 3)
#' print(result_hierarchical$cluster_information)
#'
#' # DBSCAN clustering
#' result_dbscan <- clustering(data = iris_data, method = "DBSCAN", eps = 0.5, minpts = 5)
#' print(result_dbscan$dbscan)
#' @export
clustering <- function(data, method = 'kmeans', metric = NULL ,
                       centers = NULL, eps = NULL, minpts = NULL) {

  # Uniform different inputs
  data <- if (is.matrix(data)) as.data.frame(data) else as.data.frame(data)
  if (!(is.matrix(data) || is.data.frame(data))) stop("Data must be a matrix or dataframe")
  # Kmean method
  if (method == 'kmeans') {
    if (is.null(centers) || centers <= 0) stop('A positive number for centers
                                               must be given.')

    kmeans_result <- kmeans(data, centers = centers)
    data$cluster <- as.factor(kmeans_result$cluster)
    pca <- prcomp(data[, -ncol(data)], scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2])
    pca_df$cluster <- data$cluster
    cluster_sizes <- as.list(table(kmeans_result$cluster))

    p <- ggplot(pca_df, aes(PC1, PC2, color = cluster, text = rownames(data))) +
      geom_point(size = 2) + theme_minimal() + labs(title = "K-means clustering
                                                    (PCA view)")

    return(list(plot = plotly::ggplotly(p), pca_result = pca_df,
                cluster_information = cluster_sizes, clusters = kmeans_result$cluster))
  }

  # Kmedoids method
  else if (method == 'kmedoids') {
    if (is.null(metric)) stop('Metric required for kmedoids.')

    kmedoids_result <- pam(data, k = centers, metric = metric, stand = FALSE)
    data$cluster <- as.factor(kmedoids_result$clustering)
    pca <- prcomp(data[, -ncol(data)], scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2])
    pca_df$cluster <- data$cluster
    #PC1 e PC2 spiegano la maggior parte della variabilità, see summary(PCA)
    p <- ggplot(pca_df, aes(PC1, PC2, color = cluster, text = rownames(data))) +
      geom_point(size = 2) + theme_minimal() + labs(title = "K-medoids clustering
                                                    (PCA view)")

    return(list(plot = plotly::ggplotly(p), pca_result = pca_df,
                cluster_information = table(kmedoids_result$clustering)))
  }
  # ATTENTION: no interactive plot for hierchical clustering method
  else if (method == 'hierarchical') {
    if (is.null(centers)) stop('A number for the centers must be given.')
    # Summary table
    d <- dist(data, method = "euclidean")
    hc <- hclust(d, method = "ward.D2")
    groups <- cutree(hc, k = centers)
    summary_table <- as.list(table(groups))
    # PCA for consistency with other methods
    data$cluster <- as.factor(groups)
    pca <- prcomp(data[, -ncol(data)], scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2])
    pca_df$cluster <- data$cluster
    # Plot PCA
    p_pca <- ggplot(pca_df, aes(PC1, PC2, color = cluster, text = rownames(data))) +
      geom_point(size = 2) + theme_minimal() + labs(title = "Hierarchical clustering (PCA view)")
    # Dendrogram plot
    dendro_plot <- factoextra::fviz_dend(hc, k = centers, rect = TRUE)

    return(list(plot = plotly::ggplotly(p_pca), dendrogram = dendro_plot, pca_result = pca_df,
      cluster_information = summary_table, clusters = groups))
  }

  # DBSCAN method
  else if (method == 'DBSCAN') {
    if (is.null(eps) || is.null(minpts)) {
      stop('Must insert eps and minpts for DBSCAN')
    }

    #DBSCAN lib call
    db_result <- dbscan::dbscan(data, eps = eps, minPts = minpts)
    data$cluster <- as.factor(db_result$cluster)
    #PCA
    pca <- prcomp(data[, -ncol(data)], scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2])
    pca_df$cluster <- data$cluster
    #plotting
    p <- ggplot(pca_df, aes(PC1, PC2, color = cluster, text = rownames(data))) +
      geom_point(size = 2) + theme_minimal() + labs(title = "DBSCAN (PCA view)")

    return(list(plot = plotly::ggplotly(p), pca_result = pca_df, dbscan = db_result))
  }

  # Invalid technique
  else {
    stop('Unsupported clustering technique')
  }
}
