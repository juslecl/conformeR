#' 1) creation of conformal groups `conf_group`.
#' 2) splitting of the original dataset using `initial_split` and balancing on `conf_group`.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom rsample initial_split training testing
#'
#' @param sce SingleCellExperiment with replicate_id, obs_condition, cell_type in colData
#' @param replicate_id column name for biological replicate (string)
#' @param obs_condition column name for observed condition (string)
#' @param cell_type column name for cell types (string)
#' @param size_train proportion for training set
#' @param size_cal proportion for calibration set of training set.
#'
#' @return list with original sce, proper training, calibration, and test sets
#' @export

data_processing <- function(sce, replicate_id, obs_condition, cell_type,
                            size_train = 0.5, size_cal = 0.25) {

  stopifnot(is(sce, "SingleCellExperiment"))
  set.seed(1234)

  if (!("matrix" %in% class(assay(sce,"logcounts")))){
    assay(sce,"logcounts") <- assay(sce,"logcounts") |> as.matrix()
  }

  colData(sce)[[replicate_id]] <- as.factor(colData(sce)[[replicate_id]])
  colData(sce)[[obs_condition]] <- as.factor(colData(sce)[[obs_condition]])
  colData(sce)[[cell_type]] <- as.factor(colData(sce)[[cell_type]])

  colData_df <- as.data.frame(colData(sce)) |>
    dplyr::mutate(
      conf_group = factor(paste0(.data[[replicate_id]], " x ", .data[[cell_type]]))
    )
  colData(sce)$conf_group <- colData_df$conf_group
  colData_df$row <- seq_len(nrow(colData_df))

#################################################################################
  #mat_out <- matrix(NA,
   #                 nrow = nrow(sce),
    #                ncol = ncol(sce),
    #                dimnames = list(rownames(sce), colnames(sce)))

  #for (g in levels(sce$conf_group)) {
  #  idx <- which(sce$conf_group == g)
   # sce_g <- sce[, idx]
    #vst_g <- varianceStabilizingTransformation(assay(sce_g,"counts"), blind = FALSE)
  #  mat_out[, idx] <- vst_g
  #}
  #assay(sce, "vst_assay") <- mat_out
#################################################################################

  split1 <- initial_split(colData_df, prop = size_train, strata = conf_group)
  remaining <- training(split1)
  test_idx <- testing(split1)$row

    split2 <- initial_split(remaining, prop = 1 - size_cal, strata = conf_group)
    proper_idx <- training(split2)$row
    cal_idx <- testing(split2)$row

    train_set <- sce[, proper_idx]
    cal_set <- sce[, cal_idx]
    test_set <- sce[, test_idx]

    return(list(
      train_set = sce_to_wide_tibble(train_set, obs_condition),
      cal_set = sce_to_wide_tibble(cal_set, obs_condition),
      test_set = sce_to_wide_tibble(test_set, obs_condition)
    ))
}





