#' Explorative analysis for preprocessed data
#'
#' This function performs different types of exploratory analysis on a preprocessed
#' gene expression matrix. Supported types include correlation heatmaps, PCA,
#' gene expression scatterplots, Q-Q plots and histograms.
#'
#' @param data (preprocessed) Matrix, ready for the analysis
#' @param analysis_type Character, c("corrHeatmap", "pca", "density", lm"(metti anche qq), "scatter")
#' @param gene1 Character, first gene for comparison
#' @param gene2 Character, second gene for comparison
#' @import ggplot2
#' @import plotly
#'
#' @returns description of the varoiuos output
#' @examples
#' mat <- matrix(rnorm(1000), nrow = 10, dimnames = list(paste0("Gene", 1:10), paste0("Sample", 1:100)))
#' explorative_analysis(mat, analysis_type = "CorrHeatmap")
#' explorative_analysis(mat, analysis_type = "pca")
#' explorative_analysis(mat, analysis_type = "scatter", gene1 = "Gene1", gene2 = "Gene2")
#' explorative_analysis(mat, analysis_type = "lm", gene1 = "Gene1", gene2 = "Gene2")
#' @export
explorative_analysis <- function(data, analysis_type = 'CorrHeatmap',
                                 gene1 = NULL, gene2 = NULL) {

  # Check if the input is a matrix or a dataFrame
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input must be a matrix or data frame.")
  }

  # Transform data into a matrix and handle not finite rows
  data <- as.matrix(data)
  data <- data[apply(data, 1, function(x) all(is.finite(x))), ]
  if (nrow(data) == 0 || ncol(data) == 0) stop("Matrix is empty after removing NA/Inf values.")

  # Correlation heatmap
  if (analysis_type == 'CorrHeatmap') {
    summary_stats <- summary(t(data))
    corr_matrix <- cor(t(data), use = 'pairwise.complete.obs') #correlation between genes
    #interactive heatmap with plotly
    p <- plotly::plot_ly(
      x = colnames(corr_matrix),
      y = rownames(corr_matrix),
      z = corr_matrix,
      type = "heatmap",
      colors = colorRamp(c("green", "white", "red")),
      colorbar = list(title = "Correlation"),
      showscale = TRUE
    ) %>%
      layout(
        title = "Gene Correlation Heatmap",
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "", autorange = "reversed")
      )
    return(list(summary = summary_stats, plot = p))
  }

  # Scatterplot
  else if (analysis_type == 'scatter') {
    if (is.null(gene1) || is.null(gene2)) stop('Provide 2 genes for the scatterplot.')
    if (!(gene1 %in% rownames(data)) || !(gene2 %in% rownames(data))) stop('Genes must
                                                                           be in matrix')

    #build a dataframe for the scatterplot
    df <- data.frame(
      expression = c(data[gene1, ], data[gene2, ]),
      gene = factor(rep(c(gene1, gene2), each = ncol(data))),
      sample = rep(colnames(data), 2)
    )

    p <- ggplot(df, aes(x = sample, y = expression, color = gene)) +
      geom_point(position = position_jitter(width = 0.1), size = 2) +
      theme_minimal() +
      labs(title = paste("Expression scatterplot:", gene1, "vs", gene2),
           x = "Sample",
           y = "Expression Level")
    return(p)
  }

  # QQ-plot and hist
  else if (analysis_type == 'lm') {
    if (is.null(gene1) || is.null(gene2)) stop('Provide 2 genes for lm analysis.')
    if (!(gene1 %in% rownames(data)) || !(gene2 %in% rownames(data))) stop('Genes must be in matrix')

    vals1 <- data[gene1, ]
    vals2 <- data[gene2, ]
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    qqnorm(vals1, main = sprintf('Q-Q Plot of %s', gene1)); qqline(vals1,
                                                                   col = 'red')
    qqnorm(vals2, main = sprintf('Q-Q Plot of %s', gene2)); qqline(vals2,
                                                                   col = 'red')
    hist(vals1, prob = TRUE, main = sprintf('Histogram of %s', gene1))
    hist(vals2, prob = TRUE, main = sprintf('Histogram of %s', gene2))
    return(invisible(NULL)) #don't return any obj
  }

  #PCA analysis
  else if (analysis_type == 'pca') {
    data_t <- t(data)
    pca_result <- prcomp(data_t, center = TRUE, scale. = TRUE)
    summary_pca <- summary(pca_result)
    pca_df <- as.data.frame(pca_result$x)

    #plot
    p <- ggplot(pca_df, aes(PC1, PC2)) +
      geom_point() +
      theme_minimal()

    relevant_genes <- head(pca_result$rotation)
    return(list(summary = summary_pca, plot = p, pca_df = pca_df, relevant_genes = relevant_genes))
  }

  # Invalid format
  else {
    stop('Unsupported type of analysis, check the manual for further information')
  }
}
