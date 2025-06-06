#' Simple viewer for clustering results
#'
#' A utility function to inspect and extract various elements from the results of
#' a clustering analysis. It allows visualization (plots, dendrograms), retrieval
#' of PCA or DBSCAN results, or summary information from clustering output.
#'
#' @param clustering_result List, returned by clustering_method()/cluster_number()
#' @param what Character, choose "plot", "pca", "summary", "dendrogram"
#' @return Console output or plot
#' @examples
#' # Example: load example clustering result (assuming you have one)
#' # clustering_result <- clustering_method(my_data)
#' # data_visualization(clustering_result, what = "plot")
#' # data_visualization(clustering_result, what = "summary")

#' @export
data_visualization <- function(clustering_result, what = 'plot') {
  # Match the argument 'what' to valid choices
  what <- match.arg(what, choices = c('plot', 'pca', 'summary',
                                      'dendrogram', 'clusters', 'dbscan'))

  if (what == 'plot') {
    # If plotly object is returned, return directly
    if (inherits(clustering_result$plot, "plotly")) {
      return(clustering_result$plot)
    }
    # If it's ggplot2 (for example), return ggplot
    else if (inherits(clustering_result$plot, "gg")) {
      return(clustering_result$plot)
    } else {
      stop("Plot is neither of class 'gg' nor 'plotly'.")
    }
  }
  else if (what == 'pca') {
    # Return the PCA results: a data frame or object with PCA scores
    return(clustering_result$pca_result)
  }
  else if (what == 'summary') {
    # Ensure the cluster information is a list
    if (!is.list(clustering_result$cluster_information)) {
      stop("'cluster_information' should be a list.")
    }
    return(clustering_result$cluster_information)
  }
  # For hierarchical clustering (dendrogram)
  else if (what == 'dendrogram') {
    return(clustering_result$dendrogram)
  }
  else if (what == "clusters") {
    # Ensure clusters is an integer vector, not NULL
    if (is.null(clustering_result$clusters)) {
      stop("'clusters' cannot be NULL.")
    }
    return(clustering_result$clusters)
  }
  # For DBSCAN clustering
  else if (what == 'dbscan') {
    return(clustering_result$dbscan)
  }
}
