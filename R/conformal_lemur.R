#' Apply conformal procedures to LEMUR neighborhoods
#'
#' Trains gene-specific classifiers using LEMUR neighborhoods and applies
#' conformal procedures to calibration and test data. Two procedures are
#' available: conformal clustering and conformal selection.
#'
#' Conformal clustering constructs prediction sets for the neighborhood
#' labels using a calibration quantile. Conformal selection computes
#' gene- and cell-specific conformal p-values and applies the
#' Benjamini-Hochberg procedure to select cells.
#'
#' @param pred_train A training \code{SingleCellExperiment} object.
#' @param pred_cal A calibration \code{SingleCellExperiment} object.
#' @param pred_test A test \code{SingleCellExperiment} object.
#' @param nei_train A data frame containing the LEMUR neighborhoods for the
#'   training data. Must contain a \code{name} column identifying the genes
#'   and a \code{neighborhood} column containing the corresponding cell sets.
#' @param nei_cal A data frame containing the LEMUR neighborhoods for the
#'   calibration data. Must contain a \code{name} column identifying the genes
#'   and a \code{neighborhood} column containing the corresponding cell sets.
#' @param what Character vector specifying which conformal procedures to run.
#'   Can contain \code{"conf_clustering"}, \code{"conf_selection"}, or both.
#' @param alpha Numeric value specifying the target miscoverage level.
#'   Also used as the significance level for the Benjamini-Hochberg procedure.
#' @param return_cq logical specifying whether score quantile
#' for conformal clustering should be returned
#' (for article simulations only)
#'
#' @return A data frame containing the conformal results. If
#'   \code{what = "conf_clustering"}, the output contains the conformal
#'   clustering prediction sets. If \code{what = "conf_selection"}, it
#'   contains the conformal selection results. If both procedures are
#'   requested, the results are joined by gene and cell.
#' @noRd

conformal_lemur <- function(pred_train, pred_cal, pred_test, nei_train,
    nei_cal, what, alpha, return_cq = FALSE) {
    genes <- nei_train$name
    softies <- stats::setNames(
        lapply(genes, \(g) train_classifier(pred_train, nei_train, g)), genes
    )
    predict_all <- \(sce) lapply(
        genes,
        \(g) predict_classifier(softies, sce, g)
    )
    pred_proba <- predict_all(pred_cal)
    pred_proba_test <- predict_all(pred_test)
    per_gene <- \(f) do.call(rbind.data.frame, lapply(seq_along(genes), f))
    n_cal <- ncol(pred_cal)
    test_cells <- colnames(pred_test)
    do_cc <- "conf_clustering" %in% what
    do_cs <- "conf_selection" %in% what

    res <- NULL
    if (do_cc) {
        res <- per_gene(\(i) .conf_clustering_gene(
            pred_proba[[i]]$pred, pred_proba_test[[i]]$pred, genes[i],
            test_cells, n_cal, alpha
        ))
    }
    if (do_cs) {
        score_cal <- .calibration_scores(
            pred_proba, nei_cal,
            colnames(pred_cal)
        )
        label_set <- per_gene(\(i) .conf_selection_gene(
            score_cal$scores_cal[score_cal$gene == genes[i]],
            pred_proba_test[[i]]$pred, genes[i], test_cells, n_cal, alpha
        ))
        res <- if (is.null(res)) {
            label_set
        } else {
            dplyr::left_join(label_set, res, by = c("gene", "cell"))
        }
    }
    .summarise_sets(res, alpha, keep_cq = do_cc && (return_cq || do_cs))
}
