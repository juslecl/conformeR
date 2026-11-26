#' Encode per-gene cell sets as a logical genes x cells matrix
#'
#' Cells in `defined_cells` are FALSE unless they belong to the gene's set;
#' all other cells are NA (results undefined for them).
#' @noRd
.sets_to_logical <- function(sets, genes, cells, defined_cells) {
    m <- matrix(NA, nrow = length(genes), ncol = length(cells),
        dimnames = list(genes, cells))
    m[, defined_cells] <- FALSE
    for (g in genes) {
        hits <- intersect(as.character(sets[[g]]), defined_cells)
        if (length(hits) > 0) m[g, hits] <- TRUE
    }
    m
}

#' rowData of the conformeR altExp: LEMUR neighborhood statistics (test set)
#' and per-gene conformal summaries
#' @noRd
.conformer_rowdata <- function(nei_test, conf, genes, assays, alpha) {
    nei <- as.data.frame(nei_test)[match(genes, nei_test$name), , drop = FALSE]
    keep <- !vapply(nei, is.list, logical(1)) & names(nei) != "name"
    stats <- nei[, keep, drop = FALSE]
    names(stats) <- paste0("lemur_", names(stats))
    rd <- S4Vectors::DataFrame(stats, row.names = genes, check.names = FALSE)
    rd$n_lemur_neighborhood <- rowSums(assays$lemur_neighborhood, na.rm = TRUE)
    if (!is.null(assays$conf_selection)) {
        rd$n_conf_selection <- rowSums(assays$conf_selection, na.rm = TRUE)
    }
    if ("conf_quantile" %in% names(conf)) {
        rd$conf_quantile <- conf$conf_quantile[match(genes, conf$gene)]
    }
    rd$cp_alpha <- rep(alpha, length(genes))
    rd
}

#' Store the LEMUR / conformal results in the input SingleCellExperiment
#'
#' @param sce the input SingleCellExperiment
#' @param out output of `lemur_wf()`
#' @param conf output of `conformal_lemur()` (one row per gene, list columns
#'   of selected cell names)
#' @param genes genes of interest
#' @param params list of analysis parameters stored in metadata
#' @noRd
.build_conformer_sce <- function(sce, out, conf, genes, params) {
    cells <- colnames(sce)
    test_cells <- colnames(out$pred_test)
    split <- stats::setNames(rep(NA_character_, length(cells)), cells)
    split[colnames(out$pred_train)] <- "train"
    split[colnames(out$pred_cal)] <- "cal"
    split[test_cells] <- "test"
    sce$conformeR_split <- factor(unname(split),
        levels = c("train", "cal", "test"))

    fit <- out$fit
    rd <- function(m) {
        rownames(m) <- colnames(fit)
        m[cells, , drop = FALSE]
    }
    SingleCellExperiment::reducedDim(sce, "LEMUR") <- rd(t(fit$embedding))
    SingleCellExperiment::reducedDim(sce, "LEMUR_UMAP") <-
        rd(SingleCellExperiment::reducedDim(fit, "fit_umap"))
    SingleCellExperiment::reducedDim(sce, "LEMUR_aligned_UMAP") <-
        rd(SingleCellExperiment::reducedDim(fit, "fit_al_umap"))

    de <- SummarizedExperiment::assay(fit, "DE")[genes, , drop = FALSE]
    colnames(de) <- colnames(fit)
    assays <- list(
        DE = as.matrix(de[, cells, drop = FALSE]),
        lemur_neighborhood = .sets_to_logical(
            stats::setNames(out$nei_test$neighborhood, out$nei_test$name),
            genes, cells, test_cells
        )
    )
    as_sets <- \(col) stats::setNames(conf[[col]], conf$gene)
    if ("inside_cs" %in% names(conf)) {
        assays$conf_selection <- .sets_to_logical(
            as_sets("inside_cs"), genes, cells, test_cells
        )
    }
    if (all(c("inside_cc", "outside_cc") %in% names(conf))) {
        assays$conf_clustering_inside <- .sets_to_logical(
            as_sets("inside_cc"), genes, cells, test_cells
        )
        assays$conf_clustering_outside <- .sets_to_logical(
            as_sets("outside_cc"), genes, cells, test_cells
        )
        size <- assays$conf_clustering_inside + assays$conf_clustering_outside
        storage.mode(size) <- "integer"
        assays$conf_set_size <- size
    }
    ae <- SingleCellExperiment::SingleCellExperiment(
        assays = assays,
        rowData = .conformer_rowdata(
            out$nei_test, conf, genes, assays, params$cp_alpha
        )
    )
    if ("conformeR" %in% SingleCellExperiment::altExpNames(sce)) {
        warning("Overwriting existing altExp \"conformeR\".", call. = FALSE)
    }
    SingleCellExperiment::altExp(sce, "conformeR") <- ae

    params$n_cells <- c(table(sce$conformeR_split))
    params$package_version <- tryCatch(
        as.character(utils::packageVersion("conformeR")),
        error = function(e) NA_character_
    )
    S4Vectors::metadata(sce)$conformeR <- params
    sce
}
