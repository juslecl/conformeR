test_that("plotter_conformal_selection returns a combined ggplot for multiple genes", {

  out <- make_mock_conformer_output(n_genes = 2, n_cells = 40, seed = 1)
  p <- plotter_conformal_selection(out, genes_to_plot = c("gene1", "gene2"))

  expect_s3_class(p, "ggplot")
})

test_that("plotter_conformal_selection works for a single gene", {

  out <- make_mock_conformer_output(n_genes = 1, n_cells = 40, seed = 2)
  p <- plotter_conformal_selection(out, genes_to_plot = "gene1")

  expect_s3_class(p, "ggplot")
})

test_that("plotter_conformal_selection handles a gene with no conformal-selected cells", {

  out <- make_mock_conformer_output(n_genes = 1, n_cells = 40, seed = 3)
  out$conf_results$inside_cs <- FALSE

  p <- plotter_conformal_selection(out, genes_to_plot = "gene1")

  expect_s3_class(p, "ggplot")
})

test_that("plotter_conformal_clustering returns a combined ggplot for multiple genes", {

  out <- make_mock_conformer_output(n_genes = 2, n_cells = 40, seed = 4)
  p <- plotter_conformal_clustering(out, genes_to_plot = c("gene1", "gene2"))

  expect_s3_class(p, "ggplot")
})

test_that("plotter_conformal_clustering works for a single gene", {

  out <- make_mock_conformer_output(n_genes = 1, n_cells = 40, seed = 5)
  p <- plotter_conformal_clustering(out, genes_to_plot = "gene1")

  expect_s3_class(p, "ggplot")
})
