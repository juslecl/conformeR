#' Stop with `msg` if `ok` is not TRUE
#' @noRd
check_arg <- function(ok, msg) {
    if (!isTRUE(ok)) stop(msg, call. = FALSE)
    invisible(TRUE)
}

#' Is `x` a single non-missing number in (lower, upper)?
#' (bounds included when `closed` is TRUE for that side)
#' @noRd
.is_number_in <- function(x, lower, upper, closed = c(FALSE, FALSE)) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x)) return(FALSE)
    above <- if (closed[1]) x >= lower else x > lower
    below <- if (closed[2]) x <= upper else x < upper
    above && below
}

#' Is `x` a single TRUE/FALSE?
#' @noRd
.is_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)

#' Is `x` a single non-empty string?
#' @noRd
.is_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

#' Validate the `what` argument
#' @noRd
.check_what <- function(what) {
    check_arg(
        is.character(what) && length(what) > 0 &&
            all(what %in% c("conf_selection", "conf_clustering")),
        "`what` must contain \"conf_selection\" and/or \"conf_clustering\"."
    )
}

#' Validate the arguments of conformeR()
#' @noRd
.check_conformeR_args <- function(sce, design_lemur, contrast_column,
    contrast_levels, n_embedding, use_assay, test_fraction_lemur,
    cp_train_frac, cp_cal_frac, cp_alpha, genes_of_interest, what, BPPARAM,
    verbose) {
    check_arg(
        methods::is(sce, "SingleCellExperiment"),
        "`sce` must be a SingleCellExperiment."
    )
    check_arg(
        !is.null(colnames(sce)) && !anyNA(colnames(sce)) &&
            !anyDuplicated(colnames(sce)),
        "`sce` must have unique, non-missing column names (cell IDs)."
    )
    check_arg(
        !is.null(rownames(sce)) && !anyNA(rownames(sce)) &&
            !anyDuplicated(rownames(sce)),
        "`sce` must have unique, non-missing row names (gene IDs)."
    )
    check_arg(
        .is_string(use_assay) &&
            use_assay %in% SummarizedExperiment::assayNames(sce),
        "`use_assay` must be the name of an assay of `sce`."
    )
    check_arg(
        inherits(design_lemur, "formula"), "`design_lemur` must be a formula."
    )
    cd_names <- colnames(SummarizedExperiment::colData(sce))
    covariates <- all.vars(design_lemur)
    check_arg(
        all(covariates %in% cd_names),
        sprintf(
            "All variables of `design_lemur` must be columns of colData(sce); missing: %s.",
            paste(setdiff(covariates, cd_names), collapse = ", ")
        )
    )
    check_arg(
        .is_string(contrast_column),
        "`contrast_column` must be a single character string."
    )
    check_arg(
        contrast_column %in% cd_names,
        sprintf("'%s' is not a column of colData(sce).", contrast_column)
    )
    check_arg(
        contrast_column %in% covariates,
        "`contrast_column` must be one of the variables of `design_lemur`."
    )
    if (!is.null(contrast_levels)) {
        vals <- unique(as.character(sce[[contrast_column]]))
        check_arg(
            is.character(contrast_levels) && length(contrast_levels) == 2L &&
                !anyNA(contrast_levels) &&
                contrast_levels[1] != contrast_levels[2] &&
                all(contrast_levels %in% vals),
            sprintf(
                "`contrast_levels` must be NULL or two distinct levels of '%s'.",
                contrast_column
            )
        )
    }
    check_arg(
        .is_number_in(n_embedding, 0, min(dim(sce)), closed = c(FALSE, FALSE)) &&
            n_embedding == round(n_embedding),
        "`n_embedding` must be a positive integer smaller than the number of genes and cells of `sce`."
    )
    check_arg(
        .is_number_in(test_fraction_lemur, 0, 1, closed = c(TRUE, FALSE)),
        "`test_fraction_lemur` must be a number in [0, 1)."
    )
    check_arg(
        .is_number_in(cp_train_frac, 0, 1),
        "`cp_train_frac` must be a number in (0, 1)."
    )
    check_arg(
        .is_number_in(cp_cal_frac, 0, 1),
        "`cp_cal_frac` must be a number in (0, 1)."
    )
    check_arg(
        .is_number_in(cp_alpha, 0, 1),
        "`cp_alpha` must be a number in (0, 1)."
    )
    check_arg(
        is.character(genes_of_interest) && length(genes_of_interest) > 0 &&
            !anyNA(genes_of_interest) &&
            all(genes_of_interest %in% rownames(sce)),
        "`genes_of_interest` must be a non-empty character vector with all genes present in `rownames(sce)`."
    )
    .check_what(what)
    check_arg(
        methods::is(BPPARAM, "BiocParallelParam"),
        "`BPPARAM` must be a BiocParallelParam object, e.g. BiocParallel::SerialParam()."
    )
    check_arg(.is_flag(verbose), "`verbose` must be TRUE or FALSE.")
    invisible(TRUE)
}

#' Warn when the calibration set is too small for conformal clustering
#' @noRd
.check_calibration_size <- function(n_cal, alpha, what) {
    if ("conf_clustering" %in% what &&
        ceiling((n_cal + 1) * (1 - alpha)) > n_cal) {
        warning(sprintf(
            paste0(
                "The calibration set (%d cells) is too small for ",
                "cp_alpha = %s: conformal clustering prediction sets will ",
                "contain both labels. Increase cp_cal_frac or cp_alpha."
            ), n_cal, format(alpha)
        ), call. = FALSE)
    }
    invisible(TRUE)
}

#' Validate a LEMUR prediction object (train/cal/test)
#' @noRd
.check_pred_sce <- function(x, arg) {
    check_arg(
        methods::is(x, "SingleCellExperiment"),
        sprintf("`%s` must be a SingleCellExperiment.", arg)
    )
    check_arg(
        ncol(x) > 0 && !is.null(colnames(x)),
        sprintf("`%s` must contain at least one named cell.", arg)
    )
    check_arg(
        "DE" %in% SummarizedExperiment::assayNames(x),
        sprintf("`%s` must contain a \"DE\" assay (see lemur::test_de()).", arg)
    )
    emb <- x$embedding
    check_arg(
        is.matrix(emb) && ncol(emb) == ncol(x),
        sprintf("`%s` must carry a LEMUR embedding (`%s$embedding`).", arg, arg)
    )
    invisible(TRUE)
}

#' Validate LEMUR neighborhoods encoded as logical vectors
#' @noRd
.check_nei <- function(nei, arg, n_cells) {
    check_arg(
        is.data.frame(nei) && all(c("name", "neighborhood") %in% names(nei)),
        sprintf("`%s` must be a data frame with columns `name` and `neighborhood`.", arg)
    )
    check_arg(
        is.list(nei$neighborhood) &&
            all(vapply(nei$neighborhood, is.logical, logical(1))) &&
            all(lengths(nei$neighborhood) == n_cells),
        sprintf(
            "`%s$neighborhood` must be a list of logical vectors with one entry per cell.",
            arg
        )
    )
    invisible(TRUE)
}

#' Validate the output of conformeR() before plotting
#' @noRd
.check_conformer_output <- function(x, genes_to_plot, assays) {
    check_arg(
        methods::is(x, "SingleCellExperiment") &&
            "conformeR" %in% SingleCellExperiment::altExpNames(x) &&
            "conformeR_split" %in% colnames(SummarizedExperiment::colData(x)) &&
            "LEMUR_aligned_UMAP" %in% SingleCellExperiment::reducedDimNames(x),
        "`conformer_output` needs to be the output of `conformeR::conformeR()`."
    )
    ae <- SingleCellExperiment::altExp(x, "conformeR")
    missing <- setdiff(assays, SummarizedExperiment::assayNames(ae))
    check_arg(
        length(missing) == 0,
        sprintf(
            "`conformer_output` does not contain %s; rerun `conformeR()` with the corresponding `what`.",
            paste(missing, collapse = ", ")
        )
    )
    check_arg(
        is.character(genes_to_plot) && length(genes_to_plot) > 0 &&
            !anyNA(genes_to_plot) && all(genes_to_plot %in% rownames(ae)),
        "`genes_to_plot` need to be in `rownames(altExp(conformer_output, \"conformeR\"))`."
    )
    check_arg(
        any(x$conformeR_split == "test", na.rm = TRUE),
        "`conformer_output` contains no test cells (e.g. after subsetting)."
    )
    invisible(TRUE)
}
