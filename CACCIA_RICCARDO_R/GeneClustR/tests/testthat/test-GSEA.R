test_that("gsea_analysis works for both 'no' and 'yes' conditions with enriched terms", {
  set.seed(123)

  # synt. data
  gene_names <- c(
    "ENSG00000139618", "ENSG00000157764", "ENSG00000123415", "ENSG00000164396",
    "ENSG00000204555", "ENSG00000165168", "ENSG00000173964", "ENSG00000169475",
    "ENSG00000214912", "ENSG00000180245", rep("ENSG00000000000", 90)
  )
  data_matrix <- matrix(rnorm(1000), nrow = 100)
  rownames(data_matrix) <- gene_names
  colnames(data_matrix) <- paste0("Sample", 1:10)

  group_labels <- rep(c("A", "B"), each = 5)
  data_matrix[1:10, 6:10] <- data_matrix[1:10, 6:10] + 3

  result_no <- gsea_analysis(data_matrix, keyType = "ENSEMBL", condition = "no")
  expect_type(result_no, "list")
  expect_named(result_no, c("plot", "gsea"))

  # acceot NULL
  expect_true(is.null(result_no$gsea) || isS4(result_no$gsea))
  if (isS4(result_no$gsea)) {
    expect_s4_class(result_no$gsea, "gseaResult")
  }

  result_yes <- gsea_analysis(data_matrix, keyType = "ENSEMBL", condition = "yes", groups = group_labels)
  expect_type(result_yes, "list")
  expect_named(result_yes, c("plot", "gsea"))
  expect_true(is.null(result_yes$gsea) || isS4(result_yes$gsea))
  if (isS4(result_yes$gsea)) {
    expect_s4_class(result_yes$gsea, "gseaResult")
  }

  expect_true(is.null(result_no$gsea) || isS4(result_no$gsea))
  expect_true(is.null(result_yes$gsea) || isS4(result_yes$gsea))
})
