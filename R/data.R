#' Example single-cell dataset
#'
#' A small SingleCellExperiment object used to demonstrate
#' the conformeR workflow.
#' Random subset from a real dataset (see source). Code for reproducing from
#' original dataset at:
#' \url{https://github.com/juslecl/conformeR-paper/tree/main/data_preparation}
#' 3k cells and 100 genes for 2 patients, 2 conditions, 3 cell-types.
#'
#' @format A \code{SingleCellExperiment} object.
#'
#' @source
#' Zhao, W., Dovas, A., Spinazzi, E., et al. (2021).
#' Deconvolution of cell type-specific drug responses in human tumor
#' tissue with single-cell RNA-seq. \emph{Genome Medicine}, 13, 82.
#' GEO accession: GSE148842.
#' @docType data
#' @keywords datasets
#' @usage data("sce_reduced")
#' @details
#' This object is a subset of the dataset from Zhao et al. (2021),
#' processed for use as an example dataset in \pkg{conformeR}.

"sce_reduced"
