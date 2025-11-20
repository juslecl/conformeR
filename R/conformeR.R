#' Conformal inference with counterfactual methods for multi-condition single-cell data DE detection
#'
#' This function computes conformal prediction intervals and FDR-adjusted
#' results for gene expression data stored in a
#' \linkS4class{SingleCellExperiment} object. Results are aggregated at `cell_type` level.
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom dplyr mutate select pull
#' @importFrom stats predict
#' @importFrom data.table rbindlist
#' @param sce A \linkS4class{SingleCellExperiment} object containing
#'   gene expression data ("logcounts" assay).
#' @param obs_condition Character scalar. The name of the column in
#'   `colData(sce)` indicating the observed condition (e.g., treatment vs control).
#' @param replicate_id Character scalar. The name of the column in
#'   `colData(sce)` identifying biological replicates (e.g. patient ID).
#' @param cell_type Character scalar. The name of the column in
#'   `colData(sce)` specifying cell types.
#' @param spacing Numeric scalar. Grid spacing for prediction intervals computation.
#'   Default is `0.01`.
#' @param size_train Numeric scalar between 0 and 1. Proportion of data used
#'   for training. Default is `0.5`.
#' @param size_cal Numeric scalar between 0 and 1. Proportion of data used
#'   for calibration. Default is `0.25`.
#' @param cores Integer. Number of parallel workers for computation.
#'   Default is `32`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{INT}{A tibble of conformal prediction intervals
#'     for each gene x cell x level of confidence.}
#'   \item{FDR}{A tibble summarizing combined FDR per gene and cell_type.}
#' }
#'
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' sce <- mock_sce_object()  # example input
#' res <- conformeR(sce, obs_condition="treatment",
#'                  replicate_id="patient_id", cell_type="celltype",
#'                  cores=4)
#' }
#'
#' @export

conformeR <- function(sce,
                       obs_condition,
                       replicate_id,
                       cell_type,
                       spacing = 0.01,
                       size_train = 0.5,
                       size_cal = .25,
                       cores = 32, cutoff = .05) {
  set.seed(123)
  param <- MulticoreParam(workers = cores, RNGseed = 123)

  # Split data
  splits <- data_processing(sce, replicate_id, obs_condition, cell_type, size_train, size_cal)
  properT0 <- subset_by_condition(splits$train_set, obs_condition, 0)
  properT1 <- subset_by_condition(splits$train_set, obs_condition, 1)
  calT0    <- subset_by_condition(splits$cal_set, obs_condition, 0)
  calT1    <- subset_by_condition(splits$cal_set, obs_condition, 1)
  test     <- splits$test_set

  groups     <- levels(splits$train_set$conf_group)
  alphas     <- seq(spacing, 1 - spacing, spacing)
  gene_names <- rownames(sce)

  # Iterate over groups
  results <- lapply(groups, function(g) {
    # Subset by group
    gsets <- list(
      T0 = properT0[properT0$conf_group == g, ],
      T1 = properT1[properT1$conf_group == g, ],
      C0 = calT0[calT0$conf_group == g, ],
      C1 = calT1[calT1$conf_group == g, ],
      Te = test[test$conf_group == g,]
    )
    stopifnot(ncol(gsets$T0) > 2, ncol(gsets$T1) > 2)

    idx0 <- which(gsets$Te[[obs_condition]] == 0)
    idx1 <- which(gsets$Te[[obs_condition]] == 1)

    # Loop over genes in parallel
    gene_pvalues <- BiocParallel::bplapply(gene_names, function(gene) {

      # Train quantile regressions
      qrT0 <- train_qr(gsets$T0, gene, gene_names)
      qrT1 <- train_qr(gsets$T1, gene, gene_names)

      # Propensity score
      ps_model <- prop_score(rbind(gsets$T0, gsets$T1), gene, gene_names, obs_condition)
      ps_cal   <- predict(ps_model, rbind(gsets$C0, gsets$C1), type = "prob") |> pull(.pred_1)
      ps_test <- predict(ps_model, test, type = "prob") |> pull(.pred_1)
      ps_test <- ifelse(ps_test==1,min(max(ps_test),1),ps_test)
      ps_test <- ifelse(ps_test==0,1/nrow(test),ps_test)
      w_cal    <- (1 - ps_cal) / ps_cal
      w_test <- (1 - ps_test) / ps_test
      wC0      <- w_cal[1:nrow(gsets$C0)]
      wC1      <- w_cal[(nrow(gsets$C0) + 1):length(w_cal)]

      # calibration scores
      scoresT0 <- compute_cqr_scores(qrT0, gsets$C0, gene, alphas)
      scoresT1 <- compute_cqr_scores(qrT1, gsets$C1, gene, alphas)

      int0 <- build_intervals(
        gsets$Te, idx0,
        qrT1, scoresT1, wC1, w_test,
        alphas, gene, gene_names, 1)

      int1 <- build_intervals(
        gsets$Te, idx1,
        qrT0, scoresT0, wC0, w_test,
        alphas, gene, gene_names, 0)

      int <- rbind(int0,int1) |>
        mutate(gene=gene) |>
        mutate(covered = sign(lower)!=sign(upper)) |>
        mutate(conf_group = g)

      tab_res <- fdr(int,g,gene,cutoff)
      tab_res
    }, BPPARAM = param)
    rbindlist(gene_pvalues)
  })
  tab_res <- rbindlist(results)
  fdr_tab <- comb_fdr(tab_res) |> select(-c(Rg,fdr))
  tab_res <- tab_res |> select(-c(covered,Rg))

  return(list(INT=tab_res,FDR=fdr_tab))
}
