#' Train a gene-specific classifier from LEMUR neighborhoods
#'
#' Trains a logistic regression classifier to predict LEMUR neighborhood
#' membership for a given gene. The classifier uses gene-level differential
#' expression and the LEMUR embedding as predictors. Principal component
#' analysis is applied to the predictors, retaining the components explaining
#' more than 90% of the total variance.
#'
#' @param pred_train A training \code{SingleCellExperiment} object containing
#' the differential-expression matrix and LEMUR embedding.
#' @param nei_train A data frame containing the LEMUR neighborhoods. Must
#' contain a \code{name} column identifying genes and a \code{neighborhood}
#' column containing the corresponding cell sets.
#' @param gene_name Character string specifying the gene for which the
#' classifier is trained.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{fit}}{The fitted binomial logistic regression model.}
#' \item{\code{pca}}{The PCA object fitted to the predictors.}
#' \item{\code{k}}{The number of principal components retained, explaining
#' more than 90% of the total variance.}
#' }
#'
#' @export

train_classifier <- function(pred_train, nei_train, gene_name){
  idx <- which(nei_train$name==gene_name)
  y <- as.factor(as.numeric(nei_train$neighborhood[[idx]]))
  expr_de <- SummarizedExperiment::assay(pred_train, "DE")[gene_name, ]
  embedding <- t(pred_train$embedding)
  colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
  X <- data.frame(
    expr_de = expr_de,
    embedding
  )
  pca <- prcomp(X, scale. = TRUE)
  var_expl <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  k <- which(var_expl > 0.9)[1]
  dat <- data.frame(
    y = y,
    pca$x[, 1:k]
  )
  fit <- glm(y ~ ., data = dat, family = "binomial")
  list(
    fit = fit,
    pca = pca,
    k = k
  )
}

#' Predict neighborhood membership for new cells
#'
#' Applies a trained gene-specific classifier to new cells. The differential
#' expression values and LEMUR embedding are projected onto the principal
#' components learned during training. The retained components are then used
#' to predict the probability of neighborhood membership for each cell.
#'
#' @param trained_classifier A named list of trained gene-specific classifiers,
#'   as returned by \code{train_classifier}. The list must contain an element
#'   corresponding to \code{gene_name}.
#' @param pred_sce A \code{SingleCellExperiment} object containing the
#'   differential-expression matrix and LEMUR embedding for the cells to
#'   classify.
#' @param gene_name Character string specifying the gene for which predictions
#'   are required.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{\code{pred}}{Predicted probability of neighborhood membership.}
#'     \item{\code{cell_id}}{Cell identifiers from \code{pred_sce}.}
#'     \item{\code{gene}}{The gene name.}
#'   }
#'
#' @export

predict_classifier <- function(trained_classifier, pred_sce, gene_name){
  expr_de <- SummarizedExperiment::assay(pred_sce, "DE")[gene_name, ]
  embedding <- t(pred_sce$embedding)
  colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
  X_new <- data.frame(
    expr_de = expr_de,
    embedding
  )
  pca_scores <- predict(
    trained_classifier[[gene_name]]$pca,
    newdata = X_new
  )

  dat_new <- data.frame(
    pca_scores[, seq_len(trained_classifier[[gene_name]]$k), drop = FALSE]
  )

  pred <- predict(
    trained_classifier[[gene_name]]$fit,
    newdata = dat_new,
    type = "response"
  )

  data.frame(
    pred = pred,
    cell_id = colnames(pred_sce),
    gene = gene_name
  )
}


.build_lemur_contrast <- function(fit, column, levels = NULL) {
  col_data <- SummarizedExperiment::colData(fit)
  if (!column %in% colnames(col_data)) {
    stop(sprintf("'%s' is not a column of colData(sce).", column))
  }
  if (is.null(levels)) {
    vals <- col_data[[column]]
    lv <- if (is.factor(vals)) levels(vals) else sort(unique(as.character(vals)))
    lv <- lv[!is.na(lv)]
    if (length(lv) != 2) {
      stop(sprintf(
        "Column '%s' has %d levels (%s) — specify contrast_levels = c(level_of_interest, reference_level) explicitly.",
        column, length(lv), paste(lv, collapse = ", ")
      ))
    }
    levels <- lv
  }
  make_cond_call <- function(lvl) {
    args <- rlang::list2(!!column := lvl)
    rlang::call2("cond", !!!args)
  }
  # unevaluated cond(column = levels[1]) - cond(column = levels[2])
  rlang::expr(!!make_cond_call(levels[1]) - !!make_cond_call(levels[2]))
}

