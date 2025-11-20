#' Subset a tibble by a condition
#'
#' This function filters rows of a tibble where a given column `condition`matches a `value`.
#'
#' @param tibble A tibble or data frame.
#' @param condition Name of the column used for filtering (string).
#' @param value Value to match in the condition column.
#'
#' @return A filtered tibble.
#' @export

subset_by_condition <- function(tibble, condition, value) {
  tibble[tibble[[condition]] == value, ]
}

#' Convert a SingleCellExperiment to a wide tibble
#'
#' This helper converts a SingleCellExperiment object into a wide-format tibble.
#' Meta-data columns (`conf_group` and the observed condition) are merged together with
#' the gene expression matrix (counts).
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr select all_of bind_cols
#' @importFrom SummarizedExperiment colData assay
#'
#' @param sce A SingleCellExperiment object.
#' @param obs_condition Column name (string) for the observed condition in `colData`.
#'
#' @return A tibble in wide format, with cells as rows and genes, `colData` as columns.
#' @export

sce_to_wide_tibble <- function(sce,obs_condition) {
  meta_df <- as_tibble(SummarizedExperiment::colData(sce)) %>%
    select(conf_group, all_of(obs_condition))

  expr_df <- as_tibble(t(SummarizedExperiment::assay(sce,"counts")))  # transpose: rows=cells, cols=genes
  colnames(expr_df) <- rownames(sce)

  bind_cols(meta_df, expr_df)
}

#' Convert prediction intervals to p-values
#'
#' Computes p-values from conformal prediction intervals.
#' For each `(cell_id, conf_group, gene)` combination, the p-value is computed as:
#' \deqn{(1 + \sum \text{covered}) / (1 + n)}
#'
#' @importFrom dplyr group_by summarise
#'
#' @param intervals_df A tibble with columns `cell_id`, `conf_group`, `gene`,
#'   and a logical column `covered` indicating whether the interval includes the value.
#'
#' @return A tibble with one row per `(cell_id, conf_group, gene)` and a p-value.
#' @export

interval_to_pval <- function(intervals_df){
  pval_df <-intervals_df |>
    group_by(cell_id,conf_group,gene) |>
    summarise(
      pvalue = (1 + sum(covered)) / (1 + n()),
      .groups = "drop"
    )
  return(pval_df)
}

#' Transforms conformeR p-values to empirical FDR for each gene \times cell
#'
#' Converts conformal intervals to p-values, then uses the `qvalue` package to
#' estimate local false discovery rates (lfdr).
#' A fallback approach is used when qvalue estimation fails.
#'
#' @importFrom dplyr mutate
#' @importFrom qvalue qvalue
#'
#' @param int A tibble of conformal intervals.
#' @param group A grouping variable (unused but kept for compatibility).
#' @param gene Name of the gene (unused, output grouped by gene anyway).
#' @param cutoff Numeric threshold defining significance.
#'
#' @return The input tibble with added columns `Rg`, the number of discoveries and `fdr`, the fdr.
#' @export

fdr <- function(int,group,gene,cutoff){
  # Transform intervals to p-values
  pval <- int |> interval_to_pval()

  # Transform p-values to q-values
  lfdr <- tryCatch({
    qvalue::qvalue(pval, pfdr = TRUE, lambda = seq(min(pval), max(pval), 0.05))$lfdr
  }, error = function(e) {
    tryCatch({
      pi0_manual <- 1 / (1 + length(pval))
      qvalue::qvalue(pval, pi0 = pi0_manual)$lfdr
    }, error = function(e2) NULL)
  })
  fdr_val <- ifelse(sum(pval < cutoff) == 0, 0, mean(lfdr[pval < cutoff]))
  Rg <- sum(pval <= cutoff)
  int <- int |> mutate(Rg=Rg,fdr=fdr_val)
  return(int)
}

#' Compute combined FDR across conformal groups
#'
#' Aggregates FDR values across cell types for each gene.
#' The combined FDR is computed as a weighted average of FDR values, with weights
#' given by the number of rejections `Rg`.
#'
#' @importFrom dplyr group_by slice_head select mutate if_else
#'
#' @param tab_res A tibble containing columns `gene`, `conf_group`, `Rg`, and `fdr`.
#'
#' @return A tibble with combined FDR per gene \times `cell_type` \times `conf_group`.
#' @export

comb_fdr <- function(tab_res){
  fdr_tab <- tab_res |>
    group_by(gene, conf_group) |>
    slice_head(n = 1) |>
    select(Rg, fdr, gene, conf_group) |>
    mutate(
      celltype = sub(".*x\\s*", "", conf_group)
    ) |>
    group_by(gene, celltype) |>
    mutate(
      comb_fdr = if_else(
        sum(Rg) == 0,
        mean(Rg * fdr),
        sum(Rg * fdr) / sum(Rg)
      ))
  fdr_tab
}
