# Helper: lightweight stand-in for a `lemur::LemurFit` object.
#
# Builds an object that behaves like a SingleCellExperiment but additionally
# exposes an `$embedding` accessor returning an (n_embedding x n_cells)
# numeric matrix.

setClass(
  "MockLemurFit",
  contains = "SingleCellExperiment",
  representation(embedding = "matrix")
)

setMethod("$", "MockLemurFit", function(x, name) {
  if (identical(name, "embedding")) return(x@embedding)
  callNextMethod()
})

#' Build a small mock LEMUR-like fit for unit tests.
#'
#' @param n_genes Number of genes (rows of the "DE" assay).
#' @param n_cells Number of cells (columns).
#' @param n_embedding Number of embedding dimensions.
#' @param seed RNG seed, for reproducibility.
#' @return A `MockLemurFit` object with a `"DE"` assay and an `$embedding`
#'   accessor, as used by `train_classifier()` / `predict_classifier()`.
make_mock_lemur_fit <- function(n_genes = 4, n_cells = 30, n_embedding = 3, seed = 1) {
  set.seed(seed)
  gene_names <- paste0("gene", seq_len(n_genes))
  cell_names <- paste0("cell", seq_len(n_cells))

  de_mat <- matrix(
    stats::rnorm(n_genes * n_cells),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(gene_names, cell_names)
  )
  embed_mat <- matrix(
    stats::rnorm(n_embedding * n_cells),
    nrow = n_embedding, ncol = n_cells,
    dimnames = list(paste0("dim", seq_len(n_embedding)), cell_names)
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(DE = de_mat)
  )
  methods::new("MockLemurFit", sce, embedding = embed_mat)
}

#' Build a mock neighborhood data frame, as returned by
#' `lemur::find_de_neighborhoods()`.
#'
#' @param gene_names Character vector of gene names.
#' @param cell_names Character vector of cell names (sets neighborhood length).
#' @param frac_in Fraction of cells assigned to the neighborhood for each gene.
#' @param seed RNG seed.
#' @return A data frame with a `name` column and a list-column `neighborhood`
#'   of logical vectors, one per gene, each with both `TRUE` and `FALSE`
#'   entries so downstream binomial classifiers always see two classes.
make_mock_neighborhoods <- function(gene_names, cell_names, frac_in = 0.5, seed = 2) {
  set.seed(seed)
  n_cells <- length(cell_names)
  neighborhood <- lapply(seq_along(gene_names), function(i) {
    in_nei <- rep(FALSE, n_cells)
    n_in <- max(1, min(n_cells - 1, round(frac_in * n_cells)))
    in_nei[sample.int(n_cells, n_in)] <- TRUE
    in_nei
  })
  df <- data.frame(name = gene_names, stringsAsFactors = FALSE)
  df$neighborhood <- neighborhood
  df
}

#' Build a small synthetic SingleCellExperiment for `data_processing()` tests.
#'
#' @param n_cells Number of cells.
#' @param seed RNG seed.
make_mock_sce_for_split <- function(n_cells = 120, seed = 42) {
  set.seed(seed)
  counts <- matrix(
    stats::rpois(10 * n_cells, lambda = 5),
    nrow = 10, ncol = n_cells,
    dimnames = list(paste0("gene", 1:10), paste0("cell", seq_len(n_cells)))
  )
  col_data <- S4Vectors::DataFrame(
    replicate_id = factor(rep(paste0("rep", 1:3), length.out = n_cells)),
    obs_condition = factor(rep(c("treated", "control"), length.out = n_cells)),
    cell_type = factor(rep(c("A", "B"), length.out = n_cells))
  )
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = col_data
  )
}

#' Build a mock `conformeR()` return value for the `plotter_*` tests.
#'
#' `plotter_conformal_selection()` / `plotter_conformal_clustering()` consume
#' the `list(fit_lemur, nei_lemur, conf_results)` structure returned by
#' `conformeR()`. Note that `nei_lemur$neighborhood` here holds, per gene, a
#' *character vector of cell names* (the raw `lemur::find_de_neighborhoods()`
#' format) -- unlike `make_mock_neighborhoods()` above, which produces
#' logical vectors for `nei_train`/`nei_cal` as consumed by
#' `train_classifier()`. Don't mix the two up.
#'
#' @param n_genes Number of genes.
#' @param n_cells Number of cells.
#' @param seed RNG seed.
#' @return A list with `fit_lemur` (a `SingleCellExperiment` with a `"DE"`
#'   assay and a `"fit_al_umap"` reducedDim), `nei_lemur`, and `conf_results`
#'   (containing both the clustering and selection columns).
make_mock_conformer_output <- function(n_genes = 2, n_cells = 40, seed = 1) {
  set.seed(seed)
  gene_names <- paste0("gene", seq_len(n_genes))
  cell_names <- paste0("cell", seq_len(n_cells))

  de_mat <- matrix(
    stats::rnorm(n_genes * n_cells),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(gene_names, cell_names)
  )
  fit <- SingleCellExperiment::SingleCellExperiment(
    assays = list(DE = de_mat),
    colData = S4Vectors::DataFrame(row.names = cell_names)
  )
  umap_mat <- matrix(
    stats::rnorm(n_cells * 2, sd = 2),
    nrow = n_cells, ncol = 2,
    dimnames = list(cell_names, c("UMAP1", "UMAP2"))
  )
  SingleCellExperiment::reducedDim(fit, "fit_al_umap") <- umap_mat

  neighborhood <- lapply(seq_len(n_genes), function(i) {
    sample(cell_names, size = round(0.5 * n_cells))
  })
  nei_lemur <- data.frame(name = gene_names, stringsAsFactors = FALSE)
  nei_lemur$neighborhood <- neighborhood

  conf_results <- do.call(rbind, lapply(gene_names, function(g) {
    data.frame(
      gene = g,
      cell = cell_names,
      inside_cc = sample(c(TRUE, FALSE), n_cells, replace = TRUE),
      outside_cc = sample(c(TRUE, FALSE), n_cells, replace = TRUE),
      inside_cs = sample(c(TRUE, FALSE), n_cells, replace = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  list(fit_lemur = fit, nei_lemur = nei_lemur, conf_results = conf_results)
}
