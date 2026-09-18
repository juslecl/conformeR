#' Add custom axes and coordinate system to a ggplot
#'
#' Creates a fixed coordinate system with custom arrowed axes and an optional
#' label. The default ggplot axes can be removed. The returned components can
#' be added to a ggplot object using \code{+}.
#'
#' @param label Optional character string to display next to the custom axes.
#' @param fontsize Font size of the label. Defaults to 7.
#' @param arrow_length Length of the custom axes arrows. Defaults to 10.
#' @param label_offset Offset of the label from the origin. Defaults to 1.
#' @param fix_coord Logical indicating whether to use fixed coordinates.
#' Defaults to \code{TRUE}.
#' @param remove_axes Logical indicating whether to remove the default ggplot
#' axes. Defaults to \code{TRUE}.
#' @param arrow_spec A \code{grid::arrow} specification for the custom axes.
#' Defaults to a closed arrow with arrows at both ends.
#' @param units Character string specifying the units used for the arrow and
#' label positions. Defaults to \code{"mm"}.
#' @param ... Additional arguments passed to \code{ggplot2::coord_fixed}.
#'
#' @return A list of ggplot components containing the coordinate system,
#' axis theme, custom axes, and optional label.
#'
#' @export
#'
small_axis <- function(label = NULL, fontsize = 7, arrow_length = 10, label_offset = 1, fix_coord = TRUE, remove_axes = TRUE,
                       arrow_spec = grid::arrow(ends = "both", type = "closed", angle = 20, length = unit(arrow_length / 7, units)),
                       units = "mm", ...){
  coord <- if(fix_coord){
    ggplot2::coord_fixed(clip = "off", ...)
  }else{
    NULL
  }
  axis_theme <- if(remove_axes){
    ggplot2::theme(axis.line = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank())
  }else{
    NULL
  }
  lines <- ggplot2::annotation_custom(grid::polylineGrob(x = ggplot2::unit(c(0, 0, arrow_length), units), y = ggplot2::unit(c(arrow_length, 0, 0), units),
                                                gp = grid::gpar(fill = "black"),
                                                arrow = arrow_spec))
  text <- if(! is.null(label)){
    ggplot2::annotation_custom(grid::textGrob(label = label, gp = grid::gpar(fontsize = fontsize),
                                     x = ggplot2::unit(label_offset, units), y = ggplot2::unit(label_offset, units), hjust = 0, vjust = 0))
  }else{
    NULL
  }
  list(coord, axis_theme, lines, text)
}

#' Truncate a numeric value to a specified number of significant digits
#'
#' Truncates the absolute value of a numeric input to the specified number
#' of significant digits while preserving its sign.
#'
#' @param x A numeric value.
#' @param digits Number of significant digits to retain. Defaults to 6.
#'
#' @return A numeric value truncated to the specified number of significant
#' digits.
#'
#' @export
signif_to_zero <- function(x, digits = 6){
  n_signif_digits <- digits - ceiling(log10(abs(x)))
  sign(x) * floor(abs(x) * 10^n_signif_digits) / 10^n_signif_digits
}

#' Create a diverging colour scale for differential expression
#'
#' Creates a diverging colour scale centred at zero for visualizing
#' differential-expression values. The limits are symmetric around zero,
#' with a configurable white region around the midpoint.
#'
#' @param abs_max Maximum absolute value used for the scale limits.
#' @param mid_width Width of the white region around zero.
#' Defaults to 0.1.
#' @param ... Additional arguments passed to
#' \code{ggplot2::scale_color_gradientn}.
#' @param oob Function used to handle values outside the scale limits.
#' Defaults to \code{scales::squish}.
#' @param limits Limits of the colour scale. Defaults to
#' \code{c(-abs_max, abs_max)}.
#' @param breaks Break points for the colour scale. By default, the minimum,
#' zero, and maximum values are used.
#'
#' @return A ggplot2 continuous colour scale.
#'
#' @export
scale_color_de_gradient <- function(abs_max, mid_width = 0.1, ..., oob = scales::squish, limits = c(-1, 1) * abs_max, breaks = c(-1, 0, 1) * signif_to_zero(abs_max, 1)){
  colors <- c(scales::muted("blue"), "white", "white", scales::muted("red"))
  values <- c(0, 0.5 - mid_width/2, 0.5 + mid_width/2, 1)
  ggplot2::scale_color_gradientn(oob = oob, limits = limits, breaks = breaks, colors = colors, values = values, ...)
}

