#' Clustering number
#'
#' prima funzion eper esploazione e grazie al grafico capire quanti k passare
#' seconda effettiva per calcorare i cluster con i metodi
#'
#' TODO: -Aggiungi examples nelle documentazioni
#'
#' @param data Matrix, cleaned
#' @param method Character, choose c(elbow, gapstats)
#' @import factoextra
#' @import cluster
#'
#' @return a plot to understand how maany k you should use
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

  if (method == 'elbow') {
    p <- fviz_nbclust(data, pam, method = "wss", k.max = K.max)
    return(p)
  }

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

  else {
    stop('Invalid method')
  }
}
#' Clustering method
#'
#' take in input a clustering method and gives a list of useful result
#'
#' @param data Matrix or data.frame, cleaned
#' @param method Character, method = c("kmeans", "hierarchical", "kmedoids", "DBSCAN")
#' @param centers Numerical, number of clusters
#' @param metric Character, metric for kmedoids ('euclidean','manhattan')
#' @param eps Numerical, for DBSCAN
#' @param minpts Numerical, for DBSCAN
#' @import dbscan
#' @return A list with interactive plot, PCA result and clustering info
#' @export
clustering <- function(data, method = 'kmeans', metric = NULL ,centers = NULL, eps = NULL, minpts = NULL) {

  # Uniform input
  data <- if (is.matrix(data)) as.data.frame(t(data)) else as.data.frame(data)
  if (!(is.matrix(data) || is.data.frame(data))) stop("Data must be a matrix or dataframe")
  #kmean method : ...
  if (method == 'kmeans') {
    if (is.null(centers) || centers <= 0) stop('A positive number for centers must be given.')

    kmeans_result <- kmeans(data, centers = centers)
    data$cluster <- as.factor(kmeans_result$cluster)
    pca <- prcomp(data[, -ncol(data)], scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2])
    pca_df$cluster <- data$cluster
    cluster_sizes <- as.list(table(kmeans_result$cluster))

    p <- ggplot(pca_df, aes(PC1, PC2, color = cluster, text = rownames(data))) +
      geom_point(size = 2) + theme_minimal() + labs(title = "K-means clustering (PCA view)")

    return(list(plot = plotly::ggplotly(p), pca_result = pca_df, cluster_information = cluster_sizes, clusters = kmeans_result$cluster))
  }

  else if (method == 'kmedoids') {
    if (is.null(metric)) stop('Metric required for kmedoids.')

    kmedoids_result <- pam(data, k = centers, metric = metric, stand = FALSE)
    data$cluster <- as.factor(kmedoids_result$clustering)
    pca <- prcomp(data[, -ncol(data)], scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2])
    pca_df$cluster <- data$cluster
    #PC1 e PC2 spiegano la maggior parte della variabilità, see summary(PCA)
    p <- ggplot(pca_df, aes(PC1, PC2, color = cluster, text = rownames(data))) +
      geom_point(size = 2) + theme_minimal() + labs(title = "K-medoids clustering (PCA view)")

    return(list(plot = plotly::ggplotly(p), pca_result = pca_df, cluster_information = table(kmedoids_result$clustering)))
  }
  #ATTENTION: no interactive plot here
  else if (method == 'hierarchical') {
    if (is.null(centers)) stop('A number for the centers must be given.')
    #to summary table
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
    #Dendrogram plot
    dendro_plot <- factoextra::fviz_dend(hc, k = centers, rect = TRUE)

    return(list(plot = plotly::ggplotly(p_pca), dendrogram = dendro_plot, pca_result = pca_df,
      cluster_information = summary_table, clusters = groups))
  }

  else if (method == 'DBSCAN') {
    if (is.null(eps) || is.null(minpts)) {
      stop('Must insert eps and minpts for DBSCAN')
    }

    db_result <- dbscan::dbscan(data, eps = eps, minPts = minpts)
    data$cluster <- as.factor(db_result$cluster)

    pca <- prcomp(data[, -ncol(data)], scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:2])
    pca_df$cluster <- data$cluster

    p <- ggplot(pca_df, aes(PC1, PC2, color = cluster, text = rownames(data))) +
      geom_point(size = 2) + theme_minimal() + labs(title = "DBSCAN (PCA view)")

    return(list(plot = plotly::ggplotly(p), pca_result = pca_df, dbscan = db_result))
  }

  else {
    stop('Unsupported clustering technique')
  }
}
