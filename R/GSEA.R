#' Gene set enrichment analysis
#'
#' Perform a gene set enrichment analysis of the input data
#'
#' @param data Matrix or Data Frame, with gene expression (genes in rows)
#' @param keyType Character, gene ID type (e.g. 'ENSEMBL', 'SYMBOL')
#' @param condition Character, 'yes' or 'no', indicates if groups are used
#' @param groups Vector of conditions, required if condition = 'yes'
#'
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import limma
#' @importFrom stats model.matrix
#' @return A list with GSEA result and plot (or NULL plot if empty result)
#' @export
gsea_analysis <- function(data, keyType = 'ENSEMBL', condition = 'no', groups = NULL) {
  # Check keyType
  keyType <- match.arg(keyType, choices = c('ACCNUM','ALIAS','ENSEMBL','ENSEMBLPROT','ENSEMBLTRANS','ENTREZID',
                                            'ENZYME','EVIDENCE','EVIDENCEALL','FLYBASE','FLYBASECG','FLYBASEPROT',
                                            'GENENAME','GO','GOALL','MAP','ONTOLOGY','ONTOLOGYALL','PATH','PMID',
                                            'REFSEQ','SYMBOL','UNIGENE','UNIPROT'))
  condition <- match.arg(condition, choices = c("yes", "no"))

  # Remove rows with all NA
  data <- data[rowSums(is.na(data)) < ncol(data), ]

  # Ensure finite values
  data <- data[apply(data, 1, function(x) all(is.finite(x))), ]

  # Case 1: No condition
  if (condition == 'no') {
    gene_means <- rowMeans(data, na.rm = TRUE)
    gene_list <- sort(gene_means, decreasing = TRUE)
    names(gene_list) <- rownames(data)

    gse <- clusterProfiler::gseGO(geneList = gene_list,
                                  ont = "BP",
                                  OrgDb = org.Hs.eg.db,
                                  keyType = keyType,
                                  nPerm = 1000,
                                  minGSSize = 10,
                                  maxGSSize = 500,
                                  pvalueCutoff = 0.05,
                                  verbose = TRUE)

    if (is.null(gse) || nrow(gse@result) == 0) {
      warning("No enriched terms found in GSEA.")
      return(list(plot = NULL, gsea = gse))
    }

    p <- clusterProfiler::cnetplot(gse, categorySize = "pvalue", showCategory = 5)
    return(list(plot = p, gsea = gse))
  }

  # Case 2: Condition provided
  if (condition == 'yes') {
    if (is.null(groups)) stop("Groups vector must be provided when condition = 'yes'")
    if (length(groups) != ncol(data)) stop("Length of 'groups' must match number of samples")

    design <- model.matrix(~ factor(groups))
    colnames(design) <- c("Intercept", "Group")

    fit <- limma::lmFit(data, design)
    fit <- limma::eBayes(fit)
    top_table <- limma::topTable(fit, coef = "Group", number = Inf, sort.by = "t")

    gene_list <- sort(top_table$t, decreasing = TRUE)
    names(gene_list) <- rownames(top_table)

    gse <- clusterProfiler::gseGO(geneList = gene_list,
                                  ont = "BP",
                                  OrgDb = org.Hs.eg.db,
                                  keyType = keyType,
                                  minGSSize = 10,
                                  maxGSSize = 500,
                                  pvalueCutoff = 0.05,
                                  verbose = TRUE)

    if (is.null(gse) || nrow(gse@result) == 0) {
      warning("No enriched terms found in GSEA.")
      return(list(plot = NULL, gsea = gse))
    }

    p <- clusterProfiler::cnetplot(gse, categorySize = "pvalue", showCategory = 5)
    return(list(plot = p, gsea = gse))
  }
}

