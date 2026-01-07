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
#' Meta-data columns (`conf_group` and the observed condition) are meFged together with
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

sce_to_wide_tibble <- function(sce, obs_condition) {
  meta_df <- as_tibble(colData(sce)) |>
    select(conf_group, all_of(obs_condition))
  expr_df <- as_tibble(t(assay(sce, "logcounts")))
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
#' @return A tibble with one row per `cell_id` and a p-value.
#' @export

interval_to_pval <- function(intervals_df) {
  group_by(intervals_df, cell_id) |>
    summarise(
      pvalue = (1 + sum(covered)) / (1 + n()),
      .groups = "drop"
    )
}


#' Transforms conformeR p-values to empirical FDR for each gene \times cell
#'
#' Converts conformal intervals to p-values, then uses the `qvalue` package to
#' estimate the q-value (qval).
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
#' @return The input tibble with added columns `Fg`, the number of discoveries and `fdr`, the fdr.
#' @export

fdr <- function(int, cutoff) {
  pval <- interval_to_pval(int)$pvalue
  qval <- tryCatch({
    qvalue(pval, pfdr = TRUE, lambda = seq(min(pval), max(pval), 0.05))$qvalues
  }, error = function(e) {
    tryCatch({
      pi0_manual <- 1 / (1 + length(pval))
      qvalue(pval, pi0 = pi0_manual)$qval
    }, error = function(e2) NULL
    )
  })
  fdr_val <- ifelse(sum(pval <= cutoff) == 0, 1, mean(qval[pval < cutoff]))
  Rg <- sum(pval <= cutoff)
  res <- data.frame(Rg = Rg, fdr = fdr_val)
  res <- res %>%
    mutate(
      fdr = case_when(
        is.na(fdr) & Rg < 0.2 * length(pval) ~ 1,
        is.na(fdr) & Rg > 0.2 * length(pval) ~ 1 / length(pval),
        TRUE ~ fdr
      )
    )
  return(res)
}

#' Compute combined FDR across conformal groups
#'
#' Aggregates FDR values across cell types for each gene.
#' The combined FDR is computed as a weighted average of FDR values, with weights
#' given by the number of rejections `Fg`.
#'
#' @importFrom dplyr group_by slice_head select mutate if_else
#'
#' @param tab_res A tibble containing columns `gene`, `conf_group`, `Fg`, and `fdr`.
#'
#' @return A tibble with combined FDR per gene \times `cell_type` \times `conf_group`.
#' @export

comb_fdr <- function(tab_res) {
  tab_res |>
    mutate(celltype = sub(".*x\\s*", "", conf_group)) |>
    group_by(gene, celltype) |> summarize(comb_fdr=ifelse(sum(Rg)==0, 1, sum(fdr*Rg)/sum(Rg)))
}

