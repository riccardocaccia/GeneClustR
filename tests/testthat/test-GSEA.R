test_that("gsea_analysis works for both 'no' and 'yes' conditions with enriched terms", {
  set.seed(123)

  # Dati sintetici con 100 geni e 10 campioni
  gene_names <- paste0("ENSG00000", 1:100)
  data_matrix <- matrix(rnorm(1000), nrow = 100)
  rownames(data_matrix) <- gene_names
  colnames(data_matrix) <- paste0("Sample", 1:10)

  # Aggiungiamo un segnale finto per il gruppo 2
  group_labels <- rep(c("A", "B"), each = 5)
  data_matrix[1:10, 6:10] <- data_matrix[1:10, 6:10] + 3  # upregulated in group B

  # Test senza condizione (condition = "no")
  result_no <- gsea_analysis(data_matrix, keyType = "ENSEMBL", condition = "no")
  expect_type(result_no, "list")
  expect_true("gsea" %in% names(result_no))
  expect_s3_class(result_no$gsea, "gseaResult")

  # Test con condizione (condition = "yes")
  result_yes <- gsea_analysis(data_matrix, keyType = "ENSEMBL", condition = "yes", groups = group_labels)
  expect_type(result_yes, "list")
  expect_true("gsea" %in% names(result_yes))
  expect_s3_class(result_yes$gsea, "gseaResult")

  # Verifica che almeno uno dei due abbia risultati (non vuoti)
  expect_true(nrow(result_no$gsea@result) > 0 || nrow(result_yes$gsea@result) > 0)
})

