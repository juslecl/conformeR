test_that("conformeR() runs end-to-end on a small simulated dataset", {
  data("sce_vignette")
  set.seed(1)
    out <- withCallingHandlers(
    conformeR(
      sce_vignette,
      design_lemur = ~ patient_id+condition,
      contrast_column = "condition",
      n_embedding = 20,
      genes_of_interest = rownames(sce_vignette)[1:3],
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
  expect_true(all(rownames(sce_vignette)[1:3] %in% out$conf_results$gene))
})
