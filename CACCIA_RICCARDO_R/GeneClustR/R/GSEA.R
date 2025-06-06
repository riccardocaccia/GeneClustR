#' Gene set enrichment analysis
#'
#' Performs GSEA on a gene expression dataset using GO terms (Biological Process
#' ontology). Supports both ungrouped (average expression) and grouped
#' (limma differential expression) analyses.
#'
#' @param data Matrix or Data Frame, with gene expression (genes in rows)
#' @param keyType Character, gene ID type (e.g. 'ENSEMBL', 'SYMBOL')
#' @param condition Character, 'yes' or 'no', indicates if groups are used in the analysis
#' @param groups Vector of conditions, required if condition = 'yes'
#'
#' @importFrom clusterProfiler gseGO cnetplot
#' @import org.Hs.eg.db
#' @importFrom limma lmFit eBayes topTable normalizeQuantiles
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @return A list with GSEA result and plot (or NULL plot if empty result)
#' @examples
#' mat <- matrix(rnorm(1000), nrow = 100, dimnames = list(paste0("ENSG", 1:100), paste0("Sample", 1:10)))
#' groups <- rep(c("Control", "Treated"), each = 5)
#' result <- gsea_analysis(mat, keyType = "ENSEMBL", condition = "yes", groups = groups)
#' @export
gsea_analysis <- function(data, keyType = 'ENSEMBL', condition = 'no', groups = NULL) {
  # Accepted keytype
  keyType <- match.arg(keyType, choices = c('ACCNUM','ALIAS','ENSEMBL','ENSEMBLPROT','ENSEMBLTRANS','ENTREZID',
                                            'ENZYME','EVIDENCE','EVIDENCEALL','FLYBASE','FLYBASECG','FLYBASEPROT',
                                            'GENENAME','GO','GOALL','MAP','ONTOLOGY','ONTOLOGYALL','PATH','PMID',
                                            'REFSEQ','SYMBOL','UNIGENE','UNIPROT'))
  condition <- match.arg(condition, choices = c("yes", "no"))

  # Row cleaning
  data <- data[rowSums(is.na(data)) < ncol(data), ]
  data <- data[apply(data, 1, function(x) all(is.finite(x))), ]

  # Gene list, two approach for different conditions
  gene_list <- if (condition == 'no') {
    gene_means <- rowMeans(data, na.rm = TRUE)
    setNames(sort(gene_means, decreasing = TRUE), rownames(data))
  } else {
    #valid group required
    if (is.null(groups)) stop("Groups vector must be provided when condition = 'yes'")
    if (length(groups) != ncol(data)) stop("Length of 'groups' must match number
                                           of samples")

    #linear model with limma
    design <- model.matrix(~ factor(groups))
    colnames(design) <- c("Intercept", "Group")
    fit <- limma::eBayes(limma::lmFit(data, design))
    tt <- limma::topTable(fit, coef = "Group", number = Inf, sort.by = "t")
    setNames(sort(tt$t, decreasing = TRUE), rownames(tt))
  }

  gene_list <- gene_list[!duplicated(names(gene_list))] #remove duplicate
  score_type <- if (all(gene_list > 0)) "pos" else "std" #gseGO score

  #gse
  gse <- tryCatch({
    clusterProfiler::gseGO(
      geneList = gene_list,
      ont = "BP", #only biological process
      OrgDb = org.Hs.eg.db, #organism
      keyType = keyType,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE,
      scoreType = score_type
    )
  }, error = function(e) {
    warning("gseGO failed: ", conditionMessage(e))
    NULL
  })

  # Se fallisce o non ci sono risultati, crea un oggetto gseaResult vuoto valido
  if (is.null(gse) || inherits(gse, "try-error") || nrow(gse@result) == 0) {
    warning("No enriched terms found in GSEA.")
    gse <- methods::new("gseaResult",
                        result = data.frame(),
                        geneSets = list(),
                        geneList = gene_list,
                        params = list(),
                        organism = "UNKNOWN",
                        setType = "UNKNOWN",
                        readable = FALSE
    )
    return(list(plot = NULL, gsea = gse))
  }

  #plot for enriched category
  plot <- tryCatch({
    clusterProfiler::cnetplot(gse, categorySize = "pvalue", showCategory = 5)
  }, error = function(e) {
    warning("Plot generation failed: ", conditionMessage(e))
    NULL
  })

  return(list(plot = plot, gsea = gse))
}

