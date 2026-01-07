#' Fit a propensity score model
#'
#' Internal helper function that trains a propensity score model using all
#' observed logcounts except the gene of interest. The model is estimated
#' with a random forest classifier from the \pkg{ranger} engine.
#'
#' @importFrom parsnip rand_forest set_engine set_mode
#' @importFrom workflows workflow add_formula add_model fit
#'
#' @param proper_set A data frame containing predictors (logcounts) and the
#'   treatment indicator \code{Tr}.
#' @param gene Character string giving the name of the target gene.
#' @param gene_names Character vector of all gene names.
#'
#' @return A fitted workflow object containing the trained propensity score model.
#' @export

prop_score_rf <- function(proper_set, gene, gene_names, obs_condition) {
  model_prop_score <- rand_forest() |>
    set_engine("ranger") |>
    set_mode("classification")

  # Build formula dynamically
  response <- paste0("as.factor(", obs_condition, ")")
  predictors <- paste(gene_names[gene_names != gene], collapse = " + ")
  f <- as.formula(paste(response, "~", predictors))

  prop_scores <- workflow() |>
    add_formula(f) |>
    add_model(model_prop_score) |>
    fit(proper_set)

  return(prop_scores)
}


prop_score_ridge <- function(proper_set, gene, gene_names, obs_condition) {
  
  model_prop_score <- linear_reg(penalty = tune(), mixture = 0) |>
    set_engine("glmnet") |>
    set_mode("classification")

  response <- paste0("as.factor(", obs_condition, ")")
  predictors <- paste(gene_names[gene_names != gene], collapse = " + ")
  f <- as.formula(paste(response, "~", predictors))

  wf <- workflow() |>
    add_formula(f) |>
    add_model(model_prop_score)

#5-folds validation
  set.seed(123)
  folds <- vfold_cv(proper_set, v = 5, strata = all_of(obs_condition))

  grid <- grid_regular(penalty(), levels = 30)

  res <- tune_grid(
    wf,
    resamples = folds,
    grid = grid,
    metrics = metric_set(roc_auc, accuracy)
  )

  best_penalty <- select_best(res, metric = "roc_auc")

  final_wf <- finalize_workflow(wf, best_penalty)
  prop_scores <- fit(final_wf, proper_set)

  return(prop_scores)
}


#' Train a quantile regression model
#'
#' Internal helper function that fits a quantile regression forest for a given
#' gene using the expression of all other genes as predictors. The model is
#' estimated using the \pkg{ranger} implementation of quantile regression.
#'
#' @importFrom ranger ranger
#'
#' @param data A data frame containing predictors (logcounts) and outcomes.
#' @param gene Character string giving the name of the target gene.
#' @param gene_names Character vector of all gene names.
#'
#' @return A fitted \code{ranger} object trained in quantile regression mode.
#' @export
train_qr <- function(data, gene, gene_names) {
  ranger(
    formula = as.formula(
      paste(gene, "~", paste(gene_names[gene_names != gene], collapse = "+"))
    ),
    data = data,
    quantreg = TRUE
  )
}

