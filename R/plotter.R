#' Plot conformal selection results on the LEMUR embedding
#'
#' Visualizes gene-level differential expression and conformal selection
#' results on the aligned LEMUR UMAP embedding. For each selected gene,
#' cells are shown according to their LEMUR neighborhood membership and
#' conformal selection status (density contours).
#'
#' @param conformer_output A list returned by \code{conformeR},
#' containing the LEMUR fit, LEMUR neighborhoods, and conformal results.
#' @param genes_to_plot Character vector specifying the genes to plot.
#'
#' @return A single combined \code{ggplot}/\code{cowplot} object showing
#' differential expression and conformal selection results for all
#' requested genes on the aligned LEMUR UMAP embedding, with one shared
#' "Conformal selection" legend for the green contours below the grid of
#' gene panels (instead of one legend per gene).
#'
#' @examples
#' data("sce_reduced")
#' out_conformer <- conformeR(sce_reduced,
#'                            design_lemur = ~ patient_id + condition,
#'                            contrast_column = "condition",
#'                            n_embedding=15,
#'                            genes_of_interest = c("ENSG00000163453",
#'                            "ENSG00000163734"))
#' plotter_conformal_selection(out_conformer,
#'                              genes_to_plot=c("ENSG00000163453",
#'                              "ENSG00000163734"))
#'
#' @export

plotter_conformal_selection <- function(conformer_output, genes_to_plot) {
  prep <- .prepare_conformal_data(conformer_output, genes_to_plot,
                                  "inside_cs")
  color <- c("Conformal selection" = "forestgreen")
  alpha <- \(de, m) scales::rescale(abs(de), to = c(0.05, 1))
  contour <- function(d) {
    cd <- dplyr::filter(d, inside_cs)
    if (nrow(cd) < 3 || length(unique(cd$umap[, 1])) < 2 ||
        length(unique(cd$umap[, 2])) < 2) return(NULL)
    list(ggplot2::geom_density_2d(
      data = cd, ggplot2::aes(color = "Conformal selection"),
      linewidth = 0.5, contour_var = "ndensity",
      show.legend = FALSE),
      ggplot2::scale_color_manual(values = color, guide = "none"))
  }
  plots <- .gene_panels(prep$data, prep$umap_fit, TRUE, alpha, contour)
  legend <- .contour_legend(color, "", legend.box.margin =
                           ggplot2::margin(t = -10, r = 0, b = 0, l = 0))
  cowplot::plot_grid(
    .panel_row(plots, nrow = ceiling(length(genes_to_plot) / 4)),
    legend, ncol = 1, rel_heights = c(1, 0.05))
}

#' Plot conformal clustering results on the LEMUR embedding
#'
#' Visualizes gene-level differential expression and conformal clustering
#' results on the aligned LEMUR UMAP embedding. Cells are displayed according
#' to their LEMUR neighborhood membership and the size of the conformal
#' prediction set. The cells inside and outside the conformal neighborhood
#' prediction set are combined into a single stacked plot, annotated
#' "Inside" / "Outside", with one shared "Set size" legend for the
#' black/grey contours.
#'
#' @param conformer_output A list returned by \code{run_replication},
#' containing the LEMUR fit, LEMUR neighborhoods, and conformal results.
#' @param genes_to_plot Character vector specifying the genes to plot.
#'
#' @return A single combined \code{ggplot}/\code{cowplot} object: the
#' "Inside" panel row stacked above the "Outside" panel row, each row
#' labeled accordingly, with one shared "Set size" legend (black = set
#' size 1, grey = the other set size observed in the data) below the
#' combined grid.
#'
#' @examples
#' data("sce_reduced")
#' out_conformer <- conformeR(sce_reduced,
#'                            design_lemur = ~ patient_id + condition,
#'                            contrast_column = "condition",
#'                            n_embedding=15,
#'                            genes_of_interest = c("ENSG00000163453",
#'                            "ENSG00000163734"))
#' plotter_conformal_clustering(out_conformer,
#'                              genes_to_plot=c("ENSG00000163453",
#'                              "ENSG00000163734"))
#'
#' @export
plotter_conformal_clustering <- function(conformer_output, genes_to_plot) {
  prep <- .prepare_conformal_data(conformer_output, genes_to_plot,
                                  c("inside_cc", "outside_cc"))
  data <- dplyr::mutate(prep$data, set_size = inside_cc + outside_cc)
  sizes <- setdiff(sort(unique(stats::na.omit(data$set_size))), 1)
  merge_02 <- all(c(0, 2) %in% sizes)
  other <- if (merge_02) c(setdiff(sizes, c(0, 2)), "0/2") else
    as.character(sizes)
  colors <- c(stats::setNames(rep("grey", length(other)), other),
              "1" = "black")
  data$set_size_label <- ifelse(merge_02 & data$set_size %in% c(0, 2),
                                "0/2", as.character(data$set_size))
  alpha <- \(de, m) pmin(abs(de) / (0.3 * m), 1)
  row <- function(inside) {
    contour <- \(d) .set_size_contour(
      dplyr::filter(d, set_size != 1 | inside_cc == .env$inside), colors)
    .panel_row(.gene_panels(data, prep$umap_fit, inside, alpha, contour))
  }
  labels <- cowplot::plot_grid(.row_label("Inside"), .row_label("Outside"),
                               ncol = 1)
  rows <- cowplot::plot_grid(row(TRUE), row(FALSE), ncol = 1)
  cowplot::plot_grid(
    cowplot::plot_grid(labels, rows, ncol = 2, rel_widths = c(0.04, 1)),
    .contour_legend(colors, "Set size"),
    ncol = 1, rel_heights = c(1, 0.06))
}
