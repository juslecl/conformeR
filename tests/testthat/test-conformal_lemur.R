build_conformal_inputs <- function(n_genes = 2, seed_offset = 0) {
  train_fit <- make_mock_lemur_fit(n_genes = n_genes, n_cells = 60, n_embedding = 3, seed = 100 + seed_offset)
  cal_fit <- make_mock_lemur_fit(n_genes = n_genes, n_cells = 40, n_embedding = 3, seed = 101 + seed_offset)
  test_fit <- make_mock_lemur_fit(n_genes = n_genes, n_cells = 20, n_embedding = 3, seed = 102 + seed_offset)

  gene_names <- rownames(train_fit)
  rownames(cal_fit) <- gene_names
  rownames(test_fit) <- gene_names

  nei_train <- make_mock_neighborhoods(gene_names, colnames(train_fit), seed = 200 + seed_offset)
  nei_cal <- make_mock_neighborhoods(gene_names, colnames(cal_fit), seed = 201 + seed_offset)

  list(
    pred_train = train_fit, pred_cal = cal_fit, pred_test = test_fit,
    nei_train = nei_train, nei_cal = nei_cal, gene_names = gene_names
  )
}

test_that("conformal_lemur returns joined clustering + selection results", {
  d <- build_conformal_inputs(n_genes = 2)

  out <- conformal_lemur(
    pred_train = d$pred_train, pred_cal = d$pred_cal, pred_test = d$pred_test,
    nei_train = d$nei_train, nei_cal = d$nei_cal,
    what = c("conf_clustering", "conf_selection"),
    alpha = 0.1
  )

  expect_s3_class(out, "data.frame")
  expect_true(all(c("gene", "inside_cc", "outside_cc", "inside_cs", "cp_alpha") %in% names(out)))
  expect_equal(nrow(out), length(d$gene_names))
  expect_setequal(unique(out$gene), d$gene_names)
  expect_type(out$inside_cc, "list")
  expect_type(out$inside_cs, "list")
  expect_type(out$outside_cc, "list")
})

test_that("conformal_lemur('conf_clustering') returns only the clustering columns", {
  d <- build_conformal_inputs(n_genes = 2, seed_offset = 10)

  out <- conformal_lemur(
    pred_train = d$pred_train, pred_cal = d$pred_cal, pred_test = d$pred_test,
    nei_train = d$nei_train, nei_cal = d$nei_cal,
    what = "conf_clustering",
    alpha = 0.1
  )

  expect_true(all(c("inside_cc", "outside_cc", "gene", "cp_alpha") %in% names(out)))
  expect_false("inside_cs" %in% names(out))
  expect_equal(nrow(out), length(d$gene_names))
})

test_that("conformal_lemur('conf_selection') returns only the selection columns", {
  skip_if_not_installed("tidyr")

  d <- build_conformal_inputs(n_genes = 2, seed_offset = 20)

  out <- conformal_lemur(
    pred_train = d$pred_train, pred_cal = d$pred_cal, pred_test = d$pred_test,
    nei_train = d$nei_train, nei_cal = d$nei_cal,
    what = "conf_selection",
    alpha = 0.1
  )

  expect_true(all(c("inside_cs", "gene", "cp_alpha") %in% names(out)))
})
