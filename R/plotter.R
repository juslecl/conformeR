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
#' @export

plotter_conformal_selection <- function(conformer_output, genes_to_plot){
  # 1. Prepare the results from main conformeR function.
  fit <- conformer_output$fit_lemur
  fit_small <- fit[genes_to_plot,]
  nei <- conformer_output$nei_lemur |>
    dplyr::filter(name %in% genes_to_plot) |>
    dplyr::mutate(
      neighborhood = purrr::map(neighborhood, ~ colnames(fit) %in% .x)
    ) |> dplyr::select(name,neighborhood) |>
    tidyr::unnest(neighborhood) |>
    dplyr::mutate(cell=rep(colnames(fit),length(genes_to_plot)))
  conf_nei_df <- conformer_output$conf_results |>
    dplyr::filter(gene %in% genes_to_plot) |>
    dplyr::left_join(nei, by=c("cell","gene"="name")) |>
    dplyr::select(gene,neighborhood,cell,inside_cs)
  umap_fit <- SingleCellExperiment::reducedDim(fit_small, "fit_al_umap") |> as.data.frame()

  # 2. Prepare data for plotting.
  de_plot_data <- tibble::as_tibble(SingleCellExperiment::colData(fit_small), rownames = "cell") %>%
    dplyr::mutate(umap = umap_fit) %>%
    dplyr::mutate(de = tibble::as_tibble(t(SummarizedExperiment::assay(fit_small, "DE")))) %>%
    tidyr::unnest(de, names_sep = "-") %>%
    tidyr::pivot_longer(starts_with("de-"), names_sep = "-", values_to = "de", names_to = c(NA, "gene")) %>%
    dplyr::mutate(gene = factor(gene)) %>%
    dplyr::left_join(conf_nei_df, by = c("gene", "cell")) %>%
    dplyr::mutate(inside = ifelse(neighborhood, "in", "out")) %>%
    dplyr::mutate(inside = factor(inside, levels = c("in", "out")))  %>%
    dplyr::mutate(inside_cs = ifelse(inside_cs, "in", "out")) %>%
    dplyr::mutate(inside_cs = factor(inside_cs, levels = c("in", "out")))

  # 3. Plot
  de_plots <- de_plot_data %>%
    dplyr::group_by(gene) %>%
    dplyr::group_map(\(data, key){

      abs_max <- max(abs(quantile(data$de, c(0.95, 0.05))))
      data$gene <- key[[1]][1]

      contour_data <- de_plot_data %>%
        dplyr::mutate(gene = as.character(gene)) %>%
        dplyr::filter(
          gene == as.character(key[[1]]),
          inside_cs == "in"
        )

      p <- ggplot2::ggplot(data, ggplot2::aes(x = umap[,1], y = umap[,2])) +
        ggplot2::geom_point(
          data = dplyr::filter(data, neighborhood),
          ggplot2::aes(
            color = de,
            alpha = scales::rescale(abs(de), to = c(0.05, 1))
          ),
          size = 0.5
        ) +
        ggplot2::scale_alpha_identity() +
        scale_color_de_gradient(abs_max, mid_width = 0.2, name = "") +
        ggnewscale::new_scale_color() +
        ggplot2::facet_grid(
          cols = ggplot2::vars(gene),
          labeller = ggplot2::labeller(
            gene = ggplot2::as_labeller(\(x) paste0(x))
          )
        ) +
        ggplot2::scale_x_continuous(limits = range(umap_fit[,1]) * 1.1) +
        ggplot2::scale_y_continuous(limits = range(umap_fit[,2]) * 1.1) +
        small_axis(
          fontsize = 18,
          xlim = range(umap_fit[,1]),
          ylim = range(umap_fit[,2]),
          arrow_length = 6
        ) +
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
          panel.background = ggplot2::element_rect(fill = "white", colour = NA),
          plot.background = ggplot2::element_rect(fill = "white", colour = NA),
          strip.background = ggplot2::element_rect(fill = "white", colour = NA),
          panel.grid = ggplot2::element_blank()
        )

      if (
        nrow(contour_data) >= 3 &&
        dplyr::n_distinct(contour_data$umap[,1]) > 1 &&
        dplyr::n_distinct(contour_data$umap[,2]) > 1
      ) {
        p <- p +
          ggplot2::geom_density_2d(
            data = contour_data,
            ggplot2::aes(color = "Conformal selection"),
            linewidth = 0.5,
            contour_var = "ndensity",
            show.legend = FALSE
          ) +
          ggplot2::scale_color_manual(
            values = c("Conformal selection" = "forestgreen"),
            guide = "none"
          )
      }
      p
    })

  de_plots_grid <- cowplot::plot_grid(
    plotlist = de_plots,
    nrow = ceiling(length(genes_to_plot)/4),
    rel_widths = rep(1.4, length(de_plots))
  )

  # Legend
  contour_legend <- cowplot::get_legend(
    ggplot2::ggplot(
      data.frame(
        x = c(1, 2),
        y = c(1, 1),
        grp = "Conformal selection"
      ),
      ggplot2::aes(x = x, y = y, color = grp, group = grp)
    ) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::scale_color_manual(
        values = c("Conformal selection" = "forestgreen"),
        name = ""
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12),
        legend.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        )
      )
  )

  cowplot::plot_grid(
    de_plots_grid,
    contour_legend,
    ncol = 1,
    rel_heights = c(1, 0.06)
  )
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
#' @export
plotter_conformal_clustering <- function(conformer_output, genes_to_plot){
  # 1. Prepare the results from main conformeR function.
  fit <- conformer_output$fit_lemur
  fit_small <- fit[genes_to_plot,]
  nei <- conformer_output$nei_lemur |>
    dplyr::filter(name %in% genes_to_plot) |>
    dplyr::mutate(
      neighborhood = purrr::map(neighborhood, ~ colnames(fit) %in% .x)
    ) |> dplyr::select(name,neighborhood) |>
    tidyr::unnest(neighborhood) |>
    dplyr::mutate(cell=rep(colnames(fit),length(genes_to_plot)))
  conf_nei_df <- conformer_output$conf_results |>
    dplyr::filter(gene %in% genes_to_plot) |>
    dplyr::left_join(nei, by=c("cell","gene"="name")) |>
    dplyr::select(gene,neighborhood,cell,inside_cc,outside_cc) |>
    dplyr::mutate(set_size = inside_cc + outside_cc)
  umap_fit <- SingleCellExperiment::reducedDim(fit_small, "fit_al_umap") |> as.data.frame()

  # Set sizes present in the data: black is always 1, grey is whatever else occurs (0 or 2).
  other_sizes <- setdiff(sort(unique(stats::na.omit(conf_nei_df$set_size))), 1)
  set_size_colors <- c(
    "1" = "black",
    stats::setNames(rep("grey", length(other_sizes)), as.character(other_sizes))
  )

  # 2. Prepare data for plotting (shared by the Inside and Outside rows).
  de_plot_data <- tibble::as_tibble(SingleCellExperiment::colData(fit_small), rownames = "cell") %>%
    dplyr::mutate(umap = umap_fit) %>%
    dplyr::mutate(de = tibble::as_tibble(t(SummarizedExperiment::assay(fit_small, "DE")))) %>%
    tidyr::unnest(de, names_sep = "-") %>%
    tidyr::pivot_longer(starts_with("de-"), names_sep = "-", values_to = "de", names_to = c(NA, "gene")) %>%
    dplyr::mutate(gene = factor(gene)) %>%
    dplyr::left_join(conf_nei_df, by = c("gene", "cell"))

  # 3. Inside row
  de_plots_in <- de_plot_data %>%
    dplyr::group_by(gene) %>%
    dplyr::group_map(\(data, key){

      abs_max <- max(abs(quantile(data$de, c(0.95, 0.05))))
      data$gene <- key[[1]][1]

      ggplot2::ggplot(data, ggplot2::aes(x = umap[,1], y = umap[,2])) +
        ggplot2::geom_point( data = data |> dplyr::filter(neighborhood),
                             ggplot2::aes(
                               color = de,
                               alpha = pmin(abs(de) / (0.3 * abs_max), 1)
                             ),
                             size = 0.5
        ) +
        ggplot2::scale_alpha_identity() +
        scale_color_de_gradient(abs_max, mid_width = 0.2, name = "") +

        ggnewscale::new_scale_color() +

        ggplot2::geom_density_2d( data = data |>
                                    dplyr::filter(set_size != 1 | inside_cc),
                                  ggplot2::aes(
                                    color = factor(set_size),
                                    group = factor(set_size)
                                  ),
                                  linewidth = 0.5,
                                  show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
          values = set_size_colors,
          name = "Set size",
          guide = "none"
        ) +

        ggplot2::facet_grid(
          cols = ggplot2::vars(gene)
        ) +

        ggplot2::scale_x_continuous(limits = range(umap_fit[,1]) * 1.1) +
        ggplot2::scale_y_continuous(limits = range(umap_fit[,2]) * 1.1) +

        small_axis(
          fontsize = 18,
          xlim = range(umap_fit[,1]),
          ylim = range(umap_fit[,2]),
          arrow_length = 6
        ) +

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
          panel.background = ggplot2::element_rect(fill = "white", colour = NA),
          plot.background = ggplot2::element_rect(fill = "white", colour = NA),
          strip.background = ggplot2::element_rect(fill = "white", colour = NA),
          panel.grid = ggplot2::element_blank()
        )}) %>%
    cowplot::plot_grid(
      plotlist = .,
      nrow = 1,
      rel_widths = rep(1.4, length(.))
    )

  # 4. Outside row
  de_plots_out <- de_plot_data %>%
    dplyr::group_by(gene) %>%
    dplyr::group_map(\(data, key){

      abs_max <- max(abs(quantile(data$de, c(0.95, 0.05))))
      data$gene <- key[[1]][1]

      ggplot2::ggplot(data, ggplot2::aes(x = umap[,1], y = umap[,2])) +
        ggplot2::geom_point( data = data |> dplyr::filter(!neighborhood),
                             ggplot2::aes(
                               color = de,
                               alpha = pmin(abs(de) / (0.3 * abs_max), 1)
                             ),
                             size = 0.5
        ) +
        ggplot2::scale_alpha_identity() +
        scale_color_de_gradient(abs_max, mid_width = 0.2, name = "") +

        ggnewscale::new_scale_color() +

        ggplot2::geom_density_2d( data = data |>
                                    dplyr::filter(set_size != 1 | !inside_cc),
                                  ggplot2::aes(
                                    color = factor(set_size),
                                    group = factor(set_size)
                                  ),
                                  linewidth = 0.5,
                                  show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
          values = set_size_colors,
          name = "Set size",
          guide = "none"
        ) +

        ggplot2::facet_grid(
          cols = ggplot2::vars(gene)
        ) +

        ggplot2::scale_x_continuous(limits = range(umap_fit[,1]) * 1.1) +
        ggplot2::scale_y_continuous(limits = range(umap_fit[,2]) * 1.1) +

        small_axis(
          fontsize = 18,
          xlim = range(umap_fit[,1]),
          ylim = range(umap_fit[,2]),
          arrow_length = 6
        ) +

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
          panel.background = ggplot2::element_rect(fill = "white", colour = NA),
          plot.background = ggplot2::element_rect(fill = "white", colour = NA),
          strip.background = ggplot2::element_rect(fill = "white", colour = NA),
          panel.grid = ggplot2::element_blank()
        )}) %>%
    cowplot::plot_grid(
      plotlist = .,
      nrow = 1,
      rel_widths = rep(1.4, length(.))
    )

  # Row labels ("Inside" / "Outside") placed to the left of each row.
  label_inside <- cowplot::ggdraw() +
    cowplot::draw_label("Inside", angle = 90, size = 14, fontface = "bold")
  label_outside <- cowplot::ggdraw() +
    cowplot::draw_label("Outside", angle = 90, size = 14, fontface = "bold")

  labels_col <- cowplot::plot_grid(
    label_inside, label_outside,
    ncol = 1
  )

  rows <- cowplot::plot_grid(
    de_plots_in, de_plots_out,
    ncol = 1
  )

  combined_plot <- cowplot::plot_grid(
    labels_col, rows,
    ncol = 2,
    rel_widths = c(0.04, 1)
  )

  # One shared "Set size" legend (black = 1, grey = the other set size found in the data).
  legend_df <- data.frame(
    x = rep(c(1, 0), times = length(set_size_colors)),
    y = rep(seq_along(set_size_colors) - 1, each = 2),
    grp = factor(rep(names(set_size_colors), each = 2), levels = names(set_size_colors))
  )

  set_size_legend <- cowplot::get_legend(
    ggplot2::ggplot(
      legend_df,
      ggplot2::aes(x = x, y = y, color = grp, group = grp)
    ) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::scale_color_manual(
        values = set_size_colors,
        name = "Set size"
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12),
        legend.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        )
      )
  )

  cowplot::plot_grid(
    combined_plot,
    set_size_legend,
    ncol = 1,
    rel_heights = c(1, 0.06)
  )
}
