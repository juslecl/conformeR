#' Add custom axes and coordinate system for the UMAPs.
#' Source:
#' https://github.com/const-ae/lemur-Paper/blob/master/notebooks/util.R
#' @noRd
.small_axis <- function(label = NULL,
                       fontsize = 7,
                       arrow_length = 10,
                       label_offset = 1,
                       fix_coord = TRUE,
                       remove_axes = TRUE,
                       units = "mm", ...) {
    coord <- if (fix_coord) {
        ggplot2::coord_fixed(clip = "off", ...)
    } else {
        NULL
    }
    axis_theme <- if (remove_axes) {
        ggplot2::theme(
            axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
        )
    } else {
        NULL
    }
    lines <- ggplot2::annotation_custom(grid::polylineGrob(
        x = ggplot2::unit(c(0, 0, arrow_length), units),
        y = ggplot2::unit(c(arrow_length, 0, 0), units),
        gp = grid::gpar(fill = "black"),
        arrow = grid::arrow(ends = "both",
                            type = "closed",
                            angle = 20,
                            length =
                              ggplot2::unit(arrow_length / 7, units))
    ))
    text <- if (!is.null(label)) {
        ggplot2::annotation_custom(grid::textGrob(
            label = label, gp = grid::gpar(fontsize = fontsize),
            x = ggplot2::unit(label_offset, units),
            y = ggplot2::unit(label_offset, units), hjust = 0, vjust = 0
        ))
    } else {
        NULL
    }
    list(coord, axis_theme, lines, text)
}

#' Truncate a numeric value to a specified number of significant digits
#' Source:
#' https://github.com/const-ae/lemur-Paper/blob/master/notebooks/util.R
#' @noRd

.signif_to_zero <- function(x, digits = 6) {
    n_signif_digits <- digits - ceiling(log10(abs(x)))
    sign(x) * floor(abs(x) * 10^n_signif_digits) / 10^n_signif_digits
}

#' Create a diverging colour scale for differential expression
#' Source:
#' https://github.com/const-ae/lemur-Paper/blob/master/notebooks/util.R
#' @noRd
.scale_color_de_gradient <- function(abs_max,
                                    mid_width = 0.1, ...,
                                    oob = scales::squish,
                                    limits = c(-1, 1) * abs_max,
                                    breaks =
                                      c(-1, 0, 1) * .signif_to_zero(abs_max, 1))
  {
    colors <- c(scales::muted("blue"), "white", "white", scales::muted("red"))
    values <- c(0, 0.5 - mid_width / 2, 0.5 + mid_width / 2, 1)
    ggplot2::scale_color_gradientn(oob = oob,
                                   limits = limits,
                                   breaks = breaks,
                                   colors = colors,
                                   values = values, ...)
}

#' Shared theme for conformal plots
#' @noRd
.conformal_plot_theme <- function() {
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 12),
    strip.text.y = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 12),
    legend.key.width = ggplot2::unit(1, "cm"),
    legend.key.height = ggplot2::unit(3, "mm"),
    panel.spacing.x = ggplot2::unit(15, "mm"),
    panel.spacing.y = ggplot2::unit(8, "mm"),
    legend.position = "bottom",
    panel.background = ggplot2::element_rect(
      fill = "white", colour = NA
    ),
    plot.background = ggplot2::element_rect(
      fill = "white", colour = NA
    ),
    strip.background = ggplot2::element_rect(
      fill = "white", colour = NA
    ),
    panel.grid = ggplot2::element_blank()
  )
}

#' Add shared UMAP axes
#' @noRd
.add_lemur_axes <- function(p, umap_fit) {
  p +
    ggplot2::scale_x_continuous(
      limits = range(umap_fit[, 1]) * 1.1
    ) +
    ggplot2::scale_y_continuous(
      limits = range(umap_fit[, 2]) * 1.1
    ) +
    .small_axis(
      fontsize = 18,
      xlim = range(umap_fit[, 1]),
      ylim = range(umap_fit[, 2]),
      arrow_length = 6
    )
}

#' Neighborhood membership to long format
#' @noRd
.membership_long <- function(df, id_col, list_cols, genes, cells) {
  df |>
    dplyr::filter(.data[[id_col]] %in% genes) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(list_cols),
                                \(col) purrr::map(col, \(x) cells %in% x))) |>
    dplyr::select(dplyr::all_of(c(id_col, list_cols))) |>
    tidyr::unnest(dplyr::all_of(list_cols)) |>
    dplyr::mutate(cell = rep(cells, length(genes)))
}

#' LEMUR-Fit, neighborhoods and conformal output get-together for plot
#' @noRd
.prepare_conformal_data <- function(conformer_output, genes, conf_cols) {
  fit <- conformer_output$fit_lemur
  fit_small <- fit[genes, ]
  cells <- colnames(fit)
  nei <- .membership_long(conformer_output$nei_lemur, "name",
                          "neighborhood", genes, cells)
  conf <- .membership_long(conformer_output$conf_results, "gene",
                           conf_cols, genes, cells) |>
    dplyr::left_join(nei, by = c("cell", "gene" = "name"))
  umap_fit <- SingleCellExperiment::reducedDim(fit_small, "fit_al_umap") |>
    as.data.frame()
  data <- tibble::as_tibble(SingleCellExperiment::colData(fit_small),
                            rownames = "cell") |>
    dplyr::mutate(umap = umap_fit, de = tibble::as_tibble(
      t(SummarizedExperiment::assay(fit_small, "DE")))) |>
    tidyr::unnest(de, names_sep = "-") |>
    tidyr::pivot_longer(tidyselect::starts_with("de-"), names_sep = "-",
                        values_to = "de", names_to = c(NA, "gene")) |>
    dplyr::mutate(gene = factor(gene)) |>
    dplyr::left_join(conf, by = c("gene", "cell"))
  list(data = data, umap_fit = umap_fit)
}

#' Onw panel per gene; `extra(d)` for additional contours
#' @noRd
.gene_panels <- function(data, umap_fit, in_nei, alpha_fun, extra) {
  dplyr::group_by(data, gene) |>
    dplyr::group_map(\(d, key) {
      abs_max <- max(abs(stats::quantile(d$de, c(0.95, 0.05))))
      d$gene <- key[[1]][1]
      p <- ggplot2::ggplot(d, ggplot2::aes(x = umap[, 1],
                                           y = umap[, 2])) +
        ggplot2::geom_point(
          data = dplyr::filter(d, neighborhood == .env$in_nei),
          ggplot2::aes(color = de, alpha = alpha_fun(de, abs_max)),
          size = 0.5) +
        ggplot2::scale_alpha_identity() +
        .scale_color_de_gradient(abs_max, mid_width = 0.2, name = "") +
        ggnewscale::new_scale_color() +
        ggplot2::facet_grid(cols = ggplot2::vars(gene)) +
        extra(d)
      .add_lemur_axes(p, umap_fit) + .conformal_plot_theme()
    })
}

#' Panels get-together
#' @noRd
.panel_row <- function(plots, nrow = 1) {
  cowplot::plot_grid(plotlist = plots, nrow = nrow,
                     rel_widths = rep(1.4, length(plots)))
}

#' Common legend for contours (`...` -> specific themes)
#' @noRd
.contour_legend <- function(colors, title, ...) {
  df <- data.frame(x = rep(c(1, 0), times = length(colors)),
                   y = rep(seq_along(colors) - 1, each = 2),
                   grp = factor(rep(names(colors), each = 2),
                                levels = names(colors)))
  cowplot::get_legend(
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = grp,
                                     group = grp)) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::scale_color_manual(values = colors, name = title) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12),
        legend.background = ggplot2::element_rect(fill = "white",
                                                  colour = NA),
        ...))
}

#' Contours according to set size (conformal clustering)
#' @noRd
.set_size_contour <- function(d, colors) {
  d$size_grp <- factor(d$set_size_label, levels = names(colors))
  list(ggplot2::geom_density_2d(
    data = d, ggplot2::aes(color = size_grp, group = size_grp),
    linewidth = 0.5, show.legend = FALSE),
    ggplot2::scale_color_manual(values = colors, name = "Set size",
                                guide = "none"))
}

#' Y-axis labels ("Inside" / "Outside") (conformal clustering)
#' @noRd
.row_label <- function(text) {
  cowplot::ggdraw() +
    cowplot::draw_label(text, angle = 90, size = 14, fontface = "bold")
}
