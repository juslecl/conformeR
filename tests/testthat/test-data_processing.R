strata_proportions <- function(sce, strat_by) {
  strata <- interaction(as.data.frame(SummarizedExperiment::colData(sce))[strat_by])
  prop.table(table(strata))
}

test_that("data_processing partitions all cells into train/cal/test with no overlap", {
  sce <- make_mock_sce_for_split()
  split <- data_processing(
    sce,
    strat_by = c("obs_condition", "cell_type"),
    size_train = 0.6, size_cal = 0.25
  )

  expect_named(split, c("train", "cal", "test"))
  expect_s4_class(split$train, "SingleCellExperiment")
  expect_s4_class(split$cal, "SingleCellExperiment")
  expect_s4_class(split$test, "SingleCellExperiment")

  n_total <- ncol(split$train) + ncol(split$cal) + ncol(split$test)
  expect_equal(n_total, ncol(sce))

  all_ids <- c(colnames(split$train), colnames(split$cal), colnames(split$test))
  expect_equal(length(all_ids), length(unique(all_ids)))
  expect_setequal(all_ids, colnames(sce))
})

test_that("data_processing errors when given a non-SingleCellExperiment object", {
  expect_error(data_processing(data.frame(x = 1:5), strat_by = "x"))
})

test_that("data_processing respects requested split proportions approximately", {
  sce <- make_mock_sce_for_split(n_cells = 300)
  split <- data_processing(sce, strat_by = "obs_condition", size_train = 0.5, size_cal = 0.2)
  n <- ncol(sce)
  expect_true(abs(ncol(split$train) / n - 0.5 * 0.8) < 0.1)
  expect_true(abs(ncol(split$cal) / n - 0.5 * 0.2) < 0.1)
  expect_true(abs(ncol(split$test) / n - 0.5) < 0.1)
})

test_that("data_processing works with a single stratification column", {
  sce <- make_mock_sce_for_split(n_cells = 90)
  split <- data_processing(sce, strat_by = "cell_type", size_train = 0.7, size_cal = 0.3)
  expect_equal(ncol(split$train) + ncol(split$cal) + ncol(split$test), ncol(sce))
})

test_that("data_processing keeps multi-column strata roughly balanced across splits", {
  sce <- make_mock_sce_for_split(n_cells = 600)
  split <- data_processing(
    sce,
    strat_by = c("obs_condition", "cell_type"),
    size_train = 0.6, size_cal = 0.25
  )

  strat_cols <- c("obs_condition", "cell_type")
  overall <- strata_proportions(sce, strat_cols)
  train_p <- strata_proportions(split$train, strat_cols)
  cal_p   <- strata_proportions(split$cal, strat_cols)
  test_p  <- strata_proportions(split$test, strat_cols)

  expect_setequal(names(train_p), names(overall))
  expect_setequal(names(cal_p), names(overall))
  expect_setequal(names(test_p), names(overall))

  tol <- 0.05
  for (lvl in names(overall)) {
    expect_true(
      abs(train_p[[lvl]] - overall[[lvl]]) < tol,
      info = sprintf("train stratum '%s': %.3f vs overall %.3f", lvl, train_p[[lvl]], overall[[lvl]])
    )
    expect_true(
      abs(cal_p[[lvl]] - overall[[lvl]]) < tol,
      info = sprintf("cal stratum '%s': %.3f vs overall %.3f", lvl, cal_p[[lvl]], overall[[lvl]])
    )
    expect_true(
      abs(test_p[[lvl]] - overall[[lvl]]) < tol,
      info = sprintf("test stratum '%s': %.3f vs overall %.3f", lvl, test_p[[lvl]], overall[[lvl]])
    )
  }
})

test_that("data_processing keeps a single stratification column balanced across splits", {
  sce <- make_mock_sce_for_split(n_cells = 300)
  split <- data_processing(sce, strat_by = "obs_condition", size_train = 0.5, size_cal = 0.2)

  overall <- strata_proportions(sce, "obs_condition")
  tol <- 0.05
  for (part in list(train = split$train, cal = split$cal, test = split$test)) {
    p <- strata_proportions(part, "obs_condition")
    expect_setequal(names(p), names(overall))
    for (lvl in names(overall)) {
      expect_true(abs(p[[lvl]] - overall[[lvl]]) < tol)
    }
  }
})
