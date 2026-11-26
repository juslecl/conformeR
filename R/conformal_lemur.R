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
#' @param pred_train A training \code{SingleCellExperiment} object (LEMUR
#'   fit subset) with a \code{"DE"} assay and a LEMUR embedding.
#' @param pred_cal A calibration \code{SingleCellExperiment} object.
#' @param pred_test A test \code{SingleCellExperiment} object.
#' @param nei_train A data frame containing the LEMUR neighborhoods for the
#'   training data. Must contain a \code{name} column identifying the genes
#'   and a \code{neighborhood} column containing, for each gene, a logical
#'   vector over the training cells.
#' @param nei_cal A data frame containing the LEMUR neighborhoods for the
#'   calibration data. Must contain a \code{name} column identifying the genes
#'   and a \code{neighborhood} column containing, for each gene, a logical
#'   vector over the calibration cells.
#' @param what Character vector specifying which conformal procedures to run.
#'   Can contain \code{"conf_clustering"}, \code{"conf_selection"}, or both.
#' @param alpha Error level in (0, 1): miscoverage level of the conformal
#'   clustering prediction sets (coverage \eqn{1 - \alpha}) and target false
#'   discovery rate of the Benjamini-Hochberg procedure used by conformal
#'   selection.
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
    .check_pred_sce(pred_train, "pred_train")
    .check_pred_sce(pred_cal, "pred_cal")
    .check_pred_sce(pred_test, "pred_test")
    .check_nei(nei_train, "nei_train", ncol(pred_train))
    .check_nei(nei_cal, "nei_cal", ncol(pred_cal))
    .check_what(what)
    check_arg(.is_number_in(alpha, 0, 1), "`alpha` must be a number in (0, 1).")
    check_arg(.is_flag(return_cq), "`return_cq` must be TRUE or FALSE.")
    genes <- nei_train$name
    check_arg(
        length(genes) > 0 && !anyDuplicated(genes) &&
            setequal(genes, nei_cal$name) &&
            all(genes %in% rownames(pred_train)) &&
            all(genes %in% rownames(pred_cal)) &&
            all(genes %in% rownames(pred_test)),
        paste0(
            "`nei_train` and `nei_cal` must describe the same unique genes, ",
            "all present in `pred_train`, `pred_cal` and `pred_test`."
        )
    )
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
