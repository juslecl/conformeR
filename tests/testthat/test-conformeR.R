test_that("conformeR() runs end-to-end on a small simulated dataset", {

  set.seed(1)
  n_genes <- 10
  n_cells <- 1000
  base_lambda <- 30
  counts <- matrix(
    stats::rpois(n_genes * n_cells, lambda = base_lambda),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(paste0("gene", seq_len(n_genes)), paste0("cell", seq_len(n_cells)))
  )
  is_treated <- rep(c(TRUE, FALSE), each = n_cells / 2)
  for (g in paste0("gene", 1:3)) {
    counts[g, is_treated] <- stats::rpois(sum(is_treated), lambda = base_lambda * 3)
  }

  col_data <- S4Vectors::DataFrame(condition = factor(ifelse(is_treated, "treated", "control")), patient_id=factor(rep(c("A","B"),n_cells/2)))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts, logcounts = log1p(counts)),
    colData = col_data
  )
  out <- withCallingHandlers(
    conformeR(
      sce,
      design_lemur = ~ patient_id+condition,
      contrast_column = "condition",
      n_embedding = 3,
      genes_of_interest = paste0("gene", 1:3),
      n_cores = 1,
      verbose = FALSE
    ),
    warning = function(w) {
      if (grepl("did not converge", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  expect_named(out, c("fit_lemur", "nei_lemur", "conf_results"))
  expect_s3_class(out$conf_results, "data.frame")
  expect_true(all(paste0("gene", 1:3) %in% out$conf_results$gene))
})
