#'Splitting of the original dataset using `initial_split` and balancing on covariates passed to `lemur`.
#'
#' @param sce SingleCellExperiment with replicate_id, obs_condition, cell_type in colData
#' @param strat_by vector of column names in colData(sce) on which we stratify.
#' @param size_train proportion for training set
#' @param size_cal proportion for calibration set of training set.
#'
#' @return list with original sce, proper training, calibration, and test sets
#' @export

data_processing <- function(sce, strat_by,
                            size_train = 0.3, size_cal = 0.25) {

  stopifnot(is(sce, "SingleCellExperiment"))
  colData_df <- as.data.frame(SummarizedExperiment::colData(sce)) |>
    dplyr::mutate(
      row = dplyr::row_number(),
      strata = interaction(!!!rlang::syms(strat_by))
    )
  split1 <- rsample::initial_split(colData_df, prop = size_train, strata = strata)
  remaining <- rsample::training(split1)
  test_idx <- rsample::testing(split1)$row

  split2 <- rsample::initial_split(remaining, prop = 1 - size_cal, strata = strata)
  proper_idx <- rsample::training(split2)$row
  cal_idx <- rsample::testing(split2)$row

  train_set <- sce[, proper_idx]
  cal_set <- sce[, cal_idx]
  test_set <- sce[, test_idx]

  return(list(
    train = train_set,
    cal = cal_set,
    test = test_set))
}






