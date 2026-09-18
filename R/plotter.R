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
#' @return A \code{ggplot} object, or a combined plot of the requested genes,
#' showing differential expression and conformal selection results on the
#' aligned LEMUR UMAP embedding.
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
    dplyr::mutate(cell=rep(colnames(fit),nrow(genes_to_plot)))
  conf_nei_df <- conformer_output$conf_results |>
    dplyr::filter(gene %in% genes_to_plot) |>
    dplyr::left_join(nei, by=c("cell","gene"="name")) |>
    dplyr::select(gene,neighborhood,cell,inside_cs)
  umap_fit <- SingleCellExperiment::reducedDim(fit_small, "fit_al_umap") |> as.data.frame()

  # 2. Prepare data for plotting.
  de_plot_data <- tibble::as_tibble(SingleCellExperiment::colData(fit_small), rownames = "cell") %>%
    dplyr::mutate(umap = umap_fit) %>%
    dplyr::mutate(de = tibble::as_tibble(t(SingleCellExperiment::assay(fit_small, "DE")))) %>%
    tidyr::unnest(de, names_sep = "-") %>%
    tidyr::pivot_longer(starts_with("de-"), names_sep = "-", values_to = "de", names_to = c(NA, "gene")) %>%
    dplyr::mutate(gene = factor(gene)) %>%
    dplyr::left_join(is_inside, by = c("gene", "cell")) %>%
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

      ggplot2::ggplot(data, ggplot2::aes(x = umap[,1], y = umap[,2])) +
        ggrastr::rasterise(
          ggplot2::geom_point(
            ggplot2::aes(
              color = de,
              alpha = scales::rescale(abs(de), to = c(0.05, 1))
            ),
            size = 0.5
          )
        ) +
        ggplot2::geom_density_2d(
          data = de_plot_data %>%
            dplyr::mutate(gene = as.character(gene)) %>%
            dplyr::filter(gene == as.character(key[[1]])) %>%
            dplyr::group_by(gene, inside) %>%
            dplyr::filter(inside_cs =="in") %>%
            dplyr::ungroup(),
          color = "forestgreen",
          linewidth = 0.5,
          contour_var = "ndensity"
        ) +
        ggplot2::scale_alpha_identity() +
        scale_color_de_gradient(abs_max, mid_width = 0.2) +
        ggplot2::facet_grid(
          ggplot2::vars(inside),
          ggplot2::vars(gene),
          labeller = ggplot2::labeller(
            inside = "",
            gene = ggplot2::as_labeller(\(x) paste0(x))
          ),
          switch = "y"
        ) +
        ggplot2::scale_x_continuous(limits = range(umap_fit[,1]) * 1.1) +
        ggplot2::scale_y_continuous(limits = range(umap_fit[,2]) * 1.1) +
        small_axis(
          fontsize = 18,
          xlim = range(umap_fit[,1]),
          ylim = range(umap_fit[,2]),
          arrow_length = 6
        ) +
        ggplot2::guides(
          color = ggplot2::guide_colorbar(
            barheight = ggplot2::unit(4, "mm"),
            barwidth = ggplot2::unit(30, "mm"),
            title = ""
          )
        ) +
        ggplot2::theme(
          strip.text.x = ggplot2::element_text(size = 12),
          strip.text.y = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 12),
          legend.key.width = ggplot2::unit(2, "cm"),
          legend.key.height = ggplot2::unit(3, "mm"),
          panel.spacing.x = ggplot2::unit(15, "mm"),
          panel.spacing.y = ggplot2::unit(8, "mm"),
          legend.position = "bottom"
        )
    }) %>%
    cowplot::plot_grid(
      plotlist = .,
      nrow = ceiling(length(genes_to_plot)/4),
      rel_widths = rep(1.4, length(.))
    )
  de_plots
}

#' Plot conformal clustering results on the LEMUR embedding
#'
#' Visualizes gene-level differential expression and conformal clustering
#' results on the aligned LEMUR UMAP embedding. Cells are displayed according
#' to their LEMUR neighborhood membership and the size of the conformal
#' prediction set. Separate plots are returned for cells inside and outside
#' the conformal neighborhood prediction set.
#'
#' @param conformer_output A list returned by \code{run_replication},
#' containing the LEMUR fit, LEMUR neighborhoods, and conformal results.
#' @param genes_to_plot Character vector specifying the genes to plot.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{plot_in}}{Combined plot for cells inside the conformal
#' neighborhood prediction set.}
#' \item{\code{plot_out}}{Combined plot for cells outside the conformal
#' neighborhood prediction set.}
#' }
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
    dplyr::mutate(cell=rep(colnames(fit),nrow(genes_to_plot)))
  conf_nei_df <- conformer_output$conf_results |>
    dplyr::filter(gene %in% genes_to_plot) |>
    dplyr::left_join(nei, by=c("cell","gene"="name")) |>
    dplyr::select(gene,neighborhood,cell,inside_cc,outside_cc)
  umap_fit <- SingleCellExperiment::reducedDim(fit_small, "fit_al_umap") |> as.data.frame()

  data_inside <- conf_nei_df |>
    dplyr::mutate(set_size=inside_cc + outside_cc) |>
    dplyr::filter((set_size==1 & inside_cc==TRUE)|set_size!=1) |>
    dplyr::mutate("Set size"=as.factor(set_size))|>
    dplyr::left_join(
      as.data.frame(umap_fit) |> rownames_to_column("cell"),
      by = "cell"
    )

  data_outside <- conf_nei_df |>
    dplyr::mutate(set_size=inside_cc + outside_cc) |>
    dplyr::filter((set_size==1 & inside_cc==FALSE)|set_size!=1) |>
    dplyr::mutate("Set size"=as.factor(set_size)) |>
    dplyr::left_join(
      as.data.frame(umap_fit) |> rownames_to_column("cell"),
      by = "cell"
    )

  de_plot_data_in <- tibble::as_tibble(SingleCellExperiment::colData(fit_small), rownames = "cell") %>%
    dplyr::mutate(umap = umap_fit) %>%
    dplyr::mutate(de = tibble::as_tibble(t(assay(fit_small, "DE")))) %>%
    tidyr::unnest(de, names_sep = "-") %>%
    tidyr::pivot_longer(starts_with("de-"), names_sep = "-", values_to = "de", names_to = c(NA, "gene")) %>%
    dplyr::mutate(gene = factor(gene)) %>%
    dplyr::left_join(data_inside, by = c("gene", "cell"))

  de_plots_in <- de_plot_data_in %>%
    dplyr::group_by(gene) %>%
    dplyr::group_map(\(data, key){

      abs_max <- max(abs(quantile(data$de, c(0.95, 0.05))))
      data$gene <- key[[1]][1]

      ggplot2::ggplot(data, ggplot2::aes(x = umap[,1], y = umap[,2])) +
        ggrastr::rasterise(
          ggplot2::geom_point( data = data |> dplyr::filter(neighborhood),
                               ggplot2::aes(
                        color = de,
                        alpha = pmin(abs(de) / (0.3 * abs_max), 1)
                      ),
                      size = 0.5
          )) +
        ggplot2::scale_alpha_identity() +
        scale_color_de_gradient(abs_max, mid_width = 0.2,name="") +

        ggplot2::guides(
          color = ggplot2::guide_colorbar(
            barheight = ggplot2::unit(5, "mm"),
            barwidth = ggplot2::unit(30, "mm"),
            title = ""
          )
        ) +

        ggnewscale::new_scale_color() +

        ggplot2::geom_density_2d( data=data |> dplyr::filter(!is.na(set_size)),
                                  ggplot2::aes(
                           color = factor(set_size),
                           group = factor(set_size)
                         ),
                         linewidth = 0.5,
                         show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "2" = "grey",
            "1" = "black"
          ),
          name = "Set size"
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
          strip.text = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 12),
          legend.key.width = ggplot2::unit(2, "cm"),
          legend.key.height  = ggplot2::unit(3, "mm"),
          panel.spacing.x = ggplot2::unit(15, "mm"),
          panel.spacing.y = ggplot2::unit(8, "mm"),
          legend.position = "bottom"
        )}) %>%
    cowplot::plot_grid(
      plotlist = .,
      nrow = 1,
      rel_widths = rep(1.4, length(.))
    )

  de_plot_data_out <- tibble::as_tibble(SingleCellExperiment::colData(fit_small), rownames = "cell") %>%
    dplyr::mutate(umap = umap_fit) %>%
    dplyr::mutate(de = tibble::as_tibble(t(assay(fit_small, "DE")))) %>%
    tidyr::unnest(de, names_sep = "-") %>%
    tidyr::pivot_longer(starts_with("de-"), names_sep = "-", values_to = "de", names_to = c(NA, "gene")) %>%
    dplyr::mutate(gene = factor(gene)) %>%
    dplyr::left_join(data_outside, by = c("gene", "cell"))

  de_plots_out <- de_plot_data_out %>%
    dplyr::group_by(gene) %>%
    dplyr::group_map(\(data, key){

      abs_max <- max(abs(quantile(data$de, c(0.95, 0.05))))
      data$gene <- key[[1]][1]

      ggplot2::ggplot(data, aes(x = umap[,1], y = umap[,2])) +
        ggrastr::rasterise(
          ggplot2::geom_point( data = data |> dplyr::filter(neighborhood),
                      ggplot2::aes(
                        color = de,
                        alpha = pmin(abs(de) / (0.3 * abs_max), 1)
                      ),
                      size = 0.5
          )) +
        ggplot2::scale_alpha_identity() +
        scale_color_de_gradient(abs_max, mid_width = 0.2,name="") +

        ggplot2::guides(
          color =ggplot2::guide_colorbar(
            barheight = ggplot2::unit(5, "mm"),
            barwidth = ggplot2::unit(30, "mm"),
            title = ""
          )
        ) +

        ggnewscale::new_scale_color() +

        ggplot2::geom_density_2d( data=data |> dplyr::filter(!is.na(set_size)),
                                  ggplot2::aes(
                           color = factor(set_size),
                           group = factor(set_size)
                         ),
                         linewidth = 0.5,
                         show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "2" = "grey",
            "1" = "black"
          ),
          name = "Set size"
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
          strip.text = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 12),
          legend.key.width = ggplot2::unit(2, "cm"),
          legend.key.height  = ggplot2::unit(3, "mm"),
          panel.spacing.x = ggplot2::unit(15, "mm"),
          panel.spacing.y = ggplot2::unit(8, "mm"),
          legend.position = "bottom"
        )}) %>%
    cowplot::plot_grid(
      plotlist = .,
      nrow = 1,
      rel_widths = rep(1.4, length(.))
    )

 return(list(plot_in=de_plots_in,plot_out=de_plots_out))
}
