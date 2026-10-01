test_that("train_classifier returns a fitted binomial glm with PCA reduction", {
  mock_fit <- make_mock_lemur_fit(n_genes = 4, n_cells = 30, n_embedding = 3)
  nei <- make_mock_neighborhoods(rownames(mock_fit), colnames(mock_fit))

  result <- train_classifier(mock_fit, nei, "gene1")

  expect_type(result, "list")
  expect_named(result, c("fit", "pca", "k"))
  expect_s3_class(result$fit, "glm")
  expect_identical(stats::family(result$fit)$family, "binomial")
  expect_s3_class(result$pca, "prcomp")
  expect_true(result$k >= 1 && result$k <= ncol(result$pca$x))
})

test_that("train_classifier errors for a gene absent from the neighborhood table", {
  mock_fit <- make_mock_lemur_fit()
  nei <- make_mock_neighborhoods(rownames(mock_fit), colnames(mock_fit))
  expect_error(train_classifier(mock_fit, nei, "not_a_gene"))
})

test_that("predict_classifier returns probabilities in [0, 1] with the expected shape", {
  train_fit <- make_mock_lemur_fit(n_genes = 3, n_cells = 40, n_embedding = 3, seed = 10)
  nei_train <- make_mock_neighborhoods(rownames(train_fit), colnames(train_fit), seed = 11)
  classifiers <- stats::setNames(
    lapply(rownames(train_fit), function(g) train_classifier(train_fit, nei_train, g)),
    rownames(train_fit)
  )

  test_fit <- make_mock_lemur_fit(n_genes = 3, n_cells = 15, n_embedding = 3, seed = 20)
  # Predictions are matched gene-by-gene, so the test object must share gene
  # names with the training object.
  rownames(test_fit) <- rownames(train_fit)

  out <- predict_classifier(classifiers, test_fit, "gene1")

  expect_s3_class(out, "data.frame")
  expect_named(out, c("pred", "cell_id", "gene"))
  expect_equal(nrow(out), ncol(test_fit))
  expect_true(all(out$pred >= 0 & out$pred <= 1))
  expect_true(all(out$gene == "gene1"))
  expect_identical(out$cell_id, colnames(test_fit))
})

test_that("predict_classifier errors for a gene with no trained classifier", {
  train_fit <- make_mock_lemur_fit(n_genes = 2, n_cells = 20, seed = 30)
  nei_train <- make_mock_neighborhoods(rownames(train_fit), colnames(train_fit), seed = 31)
  classifiers <- list(gene1 = train_classifier(train_fit, nei_train, "gene1"))

  expect_error(predict_classifier(classifiers, train_fit, "gene2"))
})

test_that(".build_lemur_contrast auto-detects a two-level contrast column", {
  col_data <- S4Vectors::DataFrame(
    condition = factor(c("treated", "control", "treated", "control"))
  )
  fake_fit <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1:8, nrow = 2)),
    colData = col_data
  )

  contrast <- conformeR:::.build_lemur_contrast(fake_fit, "condition")
  txt <- paste(deparse(contrast), collapse = " ")

  expect_true(is.language(contrast))
  expect_match(txt, "treated")
  expect_match(txt, "control")
})

test_that(".build_lemur_contrast errors when the column does not exist", {
  fake_fit <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1:8, nrow = 2))
  )
  expect_error(
    conformeR:::.build_lemur_contrast(fake_fit, "not_a_column"),
    "not a column"
  )
})

test_that(".build_lemur_contrast errors with >2 levels and no explicit levels given", {
  col_data <- S4Vectors::DataFrame(condition = factor(c("a", "b", "c", "a")))
  fake_fit <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1:8, nrow = 2)),
    colData = col_data
  )
  expect_error(
    conformeR:::.build_lemur_contrast(fake_fit, "condition"),
    "levels"
  )
})

test_that(".build_lemur_contrast respects explicitly supplied levels", {
  col_data <- S4Vectors::DataFrame(condition = factor(c("a", "b", "c", "a")))
  fake_fit <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1:8, nrow = 2)),
    colData = col_data
  )
  contrast <- conformeR:::.build_lemur_contrast(fake_fit, "condition", levels = c("b", "c"))
  txt <- paste(deparse(contrast), collapse = " ")

  expect_match(txt, "\"b\"")
  expect_match(txt, "\"c\"")
})
