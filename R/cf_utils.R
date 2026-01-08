#' Compute normalized weights and conformal threshold eta(x)
#'
#' Internal helper function used in weighted split-CQR to normalize calibration
#' and test weights and compute the conformal threshold \eqn{\eta(x)}.
#'
#' @param scores Numeric matrix of conformal scores from the calibration set (output of `compute_cqr_scores`) for each level in `alphas`.
#' @param weights_cal Numeric vector of calibration weights.
#' @param weights_test Numeric vector, weight associated with the test point.
#' @param alphas Numeric vector of miscoverage levels in (0,1).
#'
#' @return A numeric matrix giving the conformal threshold \eqn{\eta(x)} with dimension `nrow(data_cal) x length(alphas)`.
#' @export
eta_conformeR <- function(scores, weights_cal, weights_test, alphas) {
  denom <- sum(weights_cal) + weights_test
  p_hat <- weights_cal %*% t(1 / denom)
  p_hat_infty <- weights_test / denom

  ord <- apply(scores, 2, order)
  score_s <- apply(scores, 2, sort)

  eta_list <- lapply(seq_along(alphas), function(a) {
    p <- p_hat[ord[, a], , drop = FALSE]
    cum_p <- apply(p, 2, cumsum)
    idx_first <- apply(cum_p, 2, function(col) which(col >= 1 - alphas[a])[1])
    eta <- ifelse(!is.na(idx_first), score_s[idx_first, a], p_hat_infty)
    eta
  })
  do.call(cbind, eta_list)
}


#' Compute conformal CQR scores
#'
#' Internal helper function that computes calibration scores for prediction
#' intervals in counterfactual CQR models. Scores are defined as the maximum
#' deviation of the observed outcome from the predicted lower and upper
#' quantiles.
#'
#' @importFrom stats predict
#' @param qr_model A quantile regression model fitted with conformalCQR.
#' @param data_cal A data frame containing calibration covariates and outcomes.
#' @param gene Character string giving the name of the target gene.
#' @param alphas Numeric vector of miscoverage levels in (0,1).
#'
#' @return A numeric matrix of conformal scores with dimension `nrow(data_cal) x length(alphas)`.
#' @export
compute_cqr_scores <- function(qr_model, data_cal, gene, alphas) { # vectorized version

  # Predict lower and upper quantiles
  q_lo  <- predict(qr_model, data = data_cal, type = "quantiles", quantiles = alphas/2)$predictions
  q_hi  <- predict(qr_model, data = data_cal, type = "quantiles", quantiles = 1 - alphas/2)$predictions # matrix of dim nrow(cal_set) x length(alphas)

  # Compute conformal scores
  scores <- pmax(q_lo - data_cal[[gene]], data_cal[[gene]] - q_hi) # returns a vector of dimension nrow(cal_set)*length(alphas)

  return(matrix(scores, nrow=nrow(data_cal), byrow=F))
}


#' Build conformal prediction intervals for difference between observed and counterfactual logcounts
#'
#' Constructs prediction intervals for the difference between observed and
#' counterfactual logcounts using weighted split-CQR. Vectorized: accepts
#' multiple test rows at once.
#'
#' @importFrom stats predict
#' @param test_data A data frame of test cells (rows = cells, cols = covariates).
#' @param qr_model A quantile regression model fitted with conformalCQR.
#' @param scores Numeric matrix of conformal scores with dimension `nrow(data_cal) x length(alphas)`.
#' @param weights_cal Numeric vector of calibration weights.
#' @param w_test Numeric vector of weights associated with the test samples.
#' @param alphas Numeric vector of miscoverage levels in (0,1).
#' @param gene Character string of the outcome gene.
#' @param gene_names Character vector of all gene names.
#' @param tr_flag Integer-vector flag (0 or 1) indicating treatment/control setting for each row of `test_data`
#'
#' @return A matrix with `nrow(test_data)*length(alphas)` with three columns ("lower","upper","cell_id","alpha").
#' @export

build_intervals <- function(test_data, idx, qr_model, scores, weights_cal,
                            w_test, alphas, gene, gene_names, tr_flag) {
  test_data <- test_data[idx, ]
  tr_flag <- rep(tr_flag, nrow(test_data))
  eta_mat <- eta_conformeR(scores, weights_cal, w_test[idx], alphas)

  q_lo <- predict(qr_model, data = test_data, type = "quantiles", quantiles = alphas/2)$predictions
  q_hi <- predict(qr_model, data = test_data, type = "quantiles", quantiles = 1 - alphas/2)$predictions

  y_obs <- test_data[[gene]]
  lower <- (1-tr_flag)*(y_obs - q_hi - eta_mat) + tr_flag*(q_lo - eta_mat - y_obs)
  upper <- (1-tr_flag)*(y_obs - q_lo + eta_mat) + tr_flag*(q_hi + eta_mat - y_obs)

  int <- as_tibble(
    cbind(
      lower = c(lower),
      upper = c(upper),
      cell_id = rep(idx, length(alphas)),
      alpha = rep(alphas, each = nrow(test_data))
    )
  )

  int_avg <- int |> group_by(alpha) |> summarize(lower_avg=mean(lower),
                                                 upper_avg=mean(upper),
                                                 alpha=alpha)
  return(int_avg)
}
