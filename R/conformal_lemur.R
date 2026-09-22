#' Apply conformal procedures to LEMUR neighborhoods
#'
#' Trains gene-specific classifiers using LEMUR neighborhoods and applies
#' conformal procedures to calibration and test data. Two procedures are
#' available: conformal clustering and conformal selection.
#'
#' Conformal clustering constructs prediction sets for the neighborhood
#' labels using a calibration quantile. Conformal selection computes
#' gene- and cell-specific conformal p-values and applies the
#' Benjamini-Hochberg procedure to select cells.
#'
#' @param pred_train A training \code{SingleCellExperiment} object.
#' @param pred_cal A calibration \code{SingleCellExperiment} object.
#' @param pred_test A test \code{SingleCellExperiment} object.
#' @param nei_train A data frame containing the LEMUR neighborhoods for the
#'   training data. Must contain a \code{name} column identifying the genes
#'   and a \code{neighborhood} column containing the corresponding cell sets.
#' @param nei_cal A data frame containing the LEMUR neighborhoods for the
#'   calibration data. Must contain a \code{name} column identifying the genes
#'   and a \code{neighborhood} column containing the corresponding cell sets.
#' @param what Character vector specifying which conformal procedures to run.
#'   Can contain \code{"conf_clustering"}, \code{"conf_selection"}, or both.
#' @param alpha Numeric value specifying the target miscoverage level.
#'   Also used as the significance level for the Benjamini-Hochberg procedure.
#' @param epsilon Numeric vector specifying the label corruption level for
#'   each gene. Used to shift the conformal p-values in the conformal
#'   selection procedure.
#'
#' @return A data frame containing the conformal results. If
#'   \code{what = "conf_clustering"}, the output contains the conformal
#'   clustering prediction sets. If \code{what = "conf_selection"}, it
#'   contains the conformal selection results. If both procedures are
#'   requested, the results are joined by gene and cell.
#'
#' @export

conformal_lemur <- function(pred_train, pred_cal, pred_test, nei_train, nei_cal, what, alpha, epsilon) {
  genes <- nei_train$name
  softies <- lapply(
    genes,
    function(gene_name) {
      train_classifier(pred_train,nei_train,gene_name)
    }
  )
  names(softies) <- genes

  pred_proba <- lapply(
    genes,
    function(gene_name) {
      predict_classifier(softies, pred_cal, gene_name)
    }
  )

  pred_proba_test <- lapply(
    genes,
    function(gene_name) {
      predict_classifier(softies, pred_test, gene_name)
    }
  )

  if ("conf_clustering" %in% what){
    n_cal <- ncol(pred_cal)

    pred_set_list <- lapply(
      seq_along(genes),
      function(row) {

        p_cal <- pred_proba[[row]]$pred

        s_cal <- -pmax(p_cal, 1 - p_cal)
        q <- sort(s_cal)[ceiling((n_cal + 1) * (1 - alpha))]

        p_test <- pred_proba_test[[row]]$pred

        cbind.data.frame(
          inside_cc  = (-p_test <= q),
          outside_cc = (-(1 - p_test) <= q)
        ) |>
          dplyr::mutate(
            cell = colnames(pred_test),
            gene = genes[row]
          )
      }
    )

    pred_set <- do.call(rbind.data.frame, pred_set_list)
  }

  if ("conf_selection" %in% what){
    pred_proba <- do.call(rbind.data.frame, pred_proba)

    gt_cal <- nei_cal |>
      tidyr::unnest(neighborhood) |>
      dplyr::mutate(cell_id=rep(colnames(pred_cal),nrow(nei_cal))) |>
      dplyr::rename("gene"="name")

    score_cal <- pred_proba |> dplyr::left_join(gt_cal, by = c("cell_id", "gene"))

    score_cal <- score_cal |> dplyr::mutate(scores_cal = 1000 *neighborhood - pred)

    scores_test <- lapply(
      seq_along(genes),
      function(row) {

        cbind.data.frame(
          scores_test = -pred_proba_test[[row]]$pred,
          gene = genes[row],
          cell = colnames(pred_test)
        )
      }
    )

    conf_pval <- lapply(
      seq_along(scores_test),
      function(gene) {
        corrup_shift <- epsilon[gene]/2
        current_gene <- unique(score_cal$gene)[gene]
        cal_scores <- score_cal$scores_cal[
          score_cal$gene == current_gene
        ]

        p <- vapply(
          scores_test[[gene]][, 1],
          function(score)
            min(1,((sum(cal_scores < score)+runif(1)*(1+sum(cal_scores == score))) /
                     (ncol(pred_cal) + 1))+corrup_shift),
          numeric(1)
        )

        cbind.data.frame(
          conformal_p_val = p,
          cell = colnames(pred_test)
        )
      }
    )

    scores_test <- do.call(rbind.data.frame, scores_test)
    label_set <- lapply(seq_along(conf_pval), function(gene) {

      df <- conf_pval[[gene]] |>
        as.data.frame() |>
        dplyr::mutate(
          inside_cs = p.adjust(
            as.numeric(conformal_p_val),
            method = "BH"
          ) < alpha,
          gene = genes[gene]
        )

      df
    }) |>
      (\(lst) do.call(rbind.data.frame, lst))()
  }

  if (("conf_clustering" %in% what) & ("conf_selection" %in% what)) return(dplyr::left_join(label_set,pred_set, by=c("gene","cell")))
  else if ("conf_clustering" %in% what) return(pred_set)
  else return(label_set)
}
