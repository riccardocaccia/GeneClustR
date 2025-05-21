#' Explorative analysis for preprocessed data
#'
#' Given a preprocessed df return an explorative matrix
#' -heatmap
#' -lm model?
#' -par(mfrow=c(2,2)) magari con qq e piu grafici
#' - scatter o violin
#'
#' @param data (preprocessed) Matrix, ready for the analysis
#' @param analysis_type Character, c("corrHeatmap", "pca", "density", lm"(metti anche qq), "scatter")
#' @param gene1 Character, primo gene per confronto
#' @param gene2 Character, gene 2 per confronto
#' @import ggplot2
#' @import plotly
#'
#' @returns description of the varoiuos output
#' @export
explorative_analysis <- function(data, analysis_type = 'CorrHeatmap', gene1 = NULL, gene2 = NULL) {

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input must be a matrix or data frame.")
  }
  data <- as.matrix(data)
  data <- data[apply(data, 1, function(x) all(is.finite(x))), ]
  if (nrow(data) == 0 || ncol(data) == 0) stop("Matrix is empty after removing NA/Inf values.")

  if (analysis_type == 'CorrHeatmap') {
    summary_stats <- summary(data)
    corr_matrix <- cor(t(data), use = 'pairwise.complete.obs')
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

  else if (analysis_type == 'scatter') {
    if (is.null(gene1) || is.null(gene2)) stop('Provide 2 genes for the scatterplot.')
    if (!(gene1 %in% rownames(data)) || !(gene2 %in% rownames(data))) stop('Genes must be in matrix')

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

  else if (analysis_type == 'lm') {
    if (is.null(gene1) || is.null(gene2)) stop('Provide 2 genes for lm analysis.')
    if (!(gene1 %in% rownames(data)) || !(gene2 %in% rownames(data))) stop('Genes must be in matrix')

    vals1 <- data[gene1, ]
    vals2 <- data[gene2, ]
    par(mar = c(4, 4, 1, 1))
    qqnorm(vals1, main = sprintf('Q-Q Plot of %s', gene1)); qqline(vals1, col = 'red')
    qqnorm(vals2, main = sprintf('Q-Q Plot of %s', gene2)); qqline(vals2, col = 'red')
    hist(vals1, prob = TRUE, main = sprintf('Histogram of %s', gene1))
    hist(vals2, prob = TRUE, main = sprintf('Histogram of %s', gene2))
    return(invisible(NULL))
  }

  else if (analysis_type == 'pca') {
    data_t <- t(data)
    pca_result <- prcomp(data_t, center = TRUE, scale. = TRUE)
    summary_pca <- summary(pca_result)
    pca_df <- as.data.frame(pca_result$x)

    p <- ggplot(pca_df, aes(PC1, PC2)) +
      geom_point() +
      theme_minimal()

    relevant_genes <- head(pca_result$rotation)
    return(list(summary = summary_pca, plot = p, pca_df = pca_df, relevant_genes = relevant_genes))
  }

  else {
    stop('Unsupported type of analysis, check the manual for further information')
  }
}
