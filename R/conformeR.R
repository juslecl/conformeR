#' LEMUR analysis for neighborhoods on predicted gene differential
#' expression with conformal methods
#'
#' Runs a LEMUR analysis, identifies differential-expression neighborhoods,
#' and applies conformal procedures to a set of genes of interest.
#' The data are split into training, calibration, and test sets based on the
#' covariates specified in the LEMUR design formula.
#'
#' The conformal procedures are applied to a subset of genes of
#' interest. Genes are processed in chunks to allow parallel computation.
#'
#' @param sce A \code{SingleCellExperiment} object containing the single-cell
#' expression data.
#' @param design_lemur A formula specifying the design used by LEMUR.
#' @param contrast_column A column name specifying the comparison used for
#' differential-expression testing in LEMUR.
#' @param contrast_levels when `contrast_column` has more than two levels,
#' specify which two are tested again each other.
#' @param n_embedding Number of dimensions in the LEMUR embedding.
#' @param use_assay Character string specifying the assay to use.
#' @param test_fraction_lemur Fraction of cells used as the LEMUR test set.
#' @param cp_train_frac Fraction of cells used for conformal training.
#' @param cp_cal_frac Fraction of cells used for conformal calibration.
#' @param cp_alpha Miscoverage level for the conformal procedures.
#' @param genes_of_interest Character vector containing the names of genes
#' to include in the conformal analysis.They need to correspond to the
#' row names of the sce.
#' @param what Character vector specifying the conformal procedure(s) to run.
#' @param n_cores Number of cores used for parallel processing.
#' @param verbose Logical indicating whether to display progress information.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{fit_lemur}}{The LEMUR test-set object restricted to the
#' genes of interest.}
#' \item{\code{nei_lemur}}{The LEMUR neighborhoods for the genes of
#' interest in the test set.}
#' \item{\code{conf_results}}{The results of the conformal procedures
#' for the genes of interest.}
#' }
#'
#' @examples
#' data("sce_reduced")
#' out_conformer <- conformeR(sce_reduced,
#'     design_lemur = ~ patient_id + condition,
#'     contrast_column = "condition",
#'     n_embedding = 15,
#'     genes_of_interest = c(
#'         "ENSG00000163453",
#'         "ENSG00000163734"
#'     )
#' )
#' out_conformer
#'
#' @export

conformeR <- function(sce, design_lemur = ~1, contrast_column,
    contrast_levels = NULL, n_embedding = 60,
    use_assay = "logcounts", test_fraction_lemur = 0.2,
    cp_train_frac = 0.5, cp_cal_frac = 0.25,
    cp_alpha = 0.05, genes_of_interest,
    what = c("conf_selection", "conf_clustering"),
    n_cores = 1, verbose = TRUE) {
    out <- lemur_wf(
        sce, design_lemur, contrast_column, contrast_levels, n_embedding,
        use_assay, test_fraction_lemur, cp_train_frac, cp_cal_frac,
        verbose = verbose
    )
    chunks <- split(
        genes_of_interest, ceiling(seq_along(genes_of_interest) / 5)
    )
    chunk_results <- BiocParallel::bplapply(
        chunks,
        \(genes) conformal_lemur(
            out$pred_train[genes, ], out$pred_cal[genes, ],
            out$pred_test[genes, ],
            .nei_as_logical(out$nei_train, genes, colnames(out$pred_train)),
            .nei_as_logical(out$nei_cal, genes, colnames(out$pred_cal)),
            what, cp_alpha
        ),
        BPPARAM = BiocParallel::MulticoreParam(workers = n_cores)
    )
    list(
        fit_lemur = out$pred_test[genes_of_interest, ],
        nei_lemur = dplyr::filter(out$nei_test, name %in% genes_of_interest),
        conf_results = dplyr::bind_rows(chunk_results)
    )
}
