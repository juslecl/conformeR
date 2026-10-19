#' Train a gene-specific classifier from LEMUR neighborhoods label
#' @noRd

train_classifier <- function(pred_train, nei_train, gene_name) {
    idx <- which(nei_train$name == gene_name)
    y <- as.factor(as.numeric(nei_train$neighborhood[[idx]]))
    expr_de <- SummarizedExperiment::assay(pred_train, "DE")[gene_name, ]
    embedding <- t(pred_train$embedding)
    colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
    X <- data.frame(
        expr_de = expr_de,
        embedding
    )
    pca <- stats::prcomp(X, scale. = TRUE)
    var_expl <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    k <- which(var_expl > 0.9)[1]
    dat <- data.frame(
        y = y,
        pca$x[, seq_len(k)]
    )
    fit <- stats::glm(y ~ ., data = dat, family = "binomial")
    list(
        fit = fit,
        pca = pca,
        k = k
    )
}

#' Predict neighborhood membership probability for new cells
#' Uses the model fitted by `train_classifier`
#' @noRd

predict_classifier <- function(trained_classifier, pred_sce, gene_name) {
    expr_de <- SummarizedExperiment::assay(pred_sce, "DE")[gene_name, ]
    embedding <- t(pred_sce$embedding)
    colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
    X_new <- data.frame(
        expr_de = expr_de,
        embedding
    )
    pca_scores <- stats::predict(
        trained_classifier[[gene_name]]$pca,
        newdata = X_new
    )

    dat_new <- data.frame(
        pca_scores[, seq_len(trained_classifier[[gene_name]]$k), drop = FALSE]
    )

    pred <- stats::predict(
        trained_classifier[[gene_name]]$fit,
        newdata = dat_new,
        type = "response"
    )

    data.frame(
        pred = pred,
        cell_id = colnames(pred_sce),
        gene = gene_name
    )
}

#' Store a UMAP of the current LEMUR embedding as a reduced dimension
#' @noRd
.add_umap <- function(fit, name) {
    SingleCellExperiment::reducedDim(fit, name) <- uwot::umap(t(fit$embedding))
    fit
}

#' Fits LEMUR standard workflow to the data.
#' Splits data for training, calibration and test sets.
#' Finds LEMUR neighborhood for each subset.
#' @noRd

lemur_wf <- function(sce, design_lemur, contrast_column, contrast_levels,
    n_embedding, use_assay, test_fraction_lemur, cp_train_frac,
    cp_cal_frac, verbose = TRUE) {
    fit <- lemur::lemur(
        sce,
        design = design_lemur, n_embedding = n_embedding,
        test_fraction = test_fraction_lemur, use_assay = use_assay,
        verbose = verbose
    ) |>
        .add_umap("fit_umap") |>
        lemur::align_harmony(verbose = verbose) |>
        .add_umap("fit_al_umap")
    contrast_expr <- .build_lemur_contrast(
        fit, contrast_column, contrast_levels
    )
    fit <- rlang::inject(lemur::test_de(fit, contrast = !!contrast_expr))
    covariates <- all.vars(design_lemur)
    preds <- data_processing(
        fit,
        strat_by = covariates, size_train = cp_train_frac,
        size_cal = cp_cal_frac
    )[c("train", "cal", "test")]
    neis <- lapply(preds, \(pred) lemur::find_de_neighborhoods(
        pred,
        group_by = glmGamPoi::vars(!!!rlang::syms(covariates)),
        test_method = "edgeR", verbose = verbose
    ))
    c(
        stats::setNames(preds, paste0("pred_", names(preds))),
        stats::setNames(neis, paste0("nei_", names(neis)))
    )
}

#' Restrict LEMUR neighborhoods to some genes and encode each neighborhood
#' as a logical vector over `cells`
#' @noRd
.nei_as_logical <- function(nei, genes, cells) {
    nei |>
        dplyr::filter(name %in% genes) |>
        dplyr::mutate(
            neighborhood = purrr::map(neighborhood, \(x) cells %in% x)
        )
}

#' Contrast parser for conformeR
#' @noRd
.build_lemur_contrast <- function(fit, column, levels = NULL) {
    col_data <- SummarizedExperiment::colData(fit)
    if (!column %in% colnames(col_data)) {
        stop(sprintf("'%s' is not a column of colData(sce).", column))
    }
    if (is.null(levels)) {
        vals <- col_data[[column]]
        lv <- if (is.factor(vals)) {
            levels(vals)
        } else {
            sort(unique(as.character(vals)))
        }
        lv <- lv[!is.na(lv)]
        if (length(lv) != 2) {
            stop(sprintf(
                "Column '%s' has %d levels (%s) - specify contrast_levels =
                c(level_of_interest, reference_level) explicitly.",
                column, length(lv), paste(lv, collapse = ", ")
            ))
        }
        levels <- lv
    }
    make_cond_call <- function(lvl) {
        args <- stats::setNames(list(lvl), column)
        rlang::call2("cond", !!!args)
    }
    rlang::expr(!!make_cond_call(levels[1]) - !!make_cond_call(levels[2]))
}


#' Conformal clustering prediction sets for one gene
#' @noRd
.conf_clustering_gene <- function(p_cal, p_test, gene, cells, n_cal, alpha) {
    s_cal <- 1 - pmax(p_cal, 1 - p_cal)
    q <- sort(s_cal)[ceiling((n_cal + 1) * (1 - alpha))]
    data.frame(
        gene = gene, cell = cells,
        inside_cc = 1 - p_test <= q, outside_cc = p_test <= q,
        conf_quantile = q
    )
}

#' Calibration scores (1000 * neighborhood - pred) for all genes
#' @noRd
.calibration_scores <- function(pred_proba, nei_cal, cal_cells) {
    gt_cal <- nei_cal |>
        tidyr::unnest(neighborhood) |>
        dplyr::mutate(cell_id = rep(cal_cells, nrow(nei_cal))) |>
        dplyr::rename("gene" = "name")
    do.call(rbind.data.frame, pred_proba) |>
        dplyr::left_join(gt_cal, by = c("cell_id", "gene")) |>
        dplyr::mutate(scores_cal = 1000 * neighborhood - pred)
}

#' Randomised conformal p-values + BH selection for one gene
#' @noRd
.conf_selection_gene <- function(cal_scores, p_test, gene, cells, n_cal,
    alpha) {
    s_test <- -p_test
    below <- vapply(s_test, \(s) sum(cal_scores < s), numeric(1))
    ties <- vapply(s_test, \(s) sum(cal_scores == s), numeric(1))
    p <- (below + stats::runif(length(s_test)) * (1 + ties)) / (n_cal + 1)
    data.frame(
        gene = gene, cell = cells,
        inside_cs = stats::p.adjust(p, method = "BH") < alpha
    )
}

#' Collapse cell-level results into per-gene cell lists
#' @noRd
.summarise_sets <- function(res, alpha, keep_cq) {
    cols <- intersect(c("inside_cs", "inside_cc", "outside_cc"), names(res))
    res <- dplyr::group_by(res, gene)
    out <- if (keep_cq) {
        dplyr::summarise(res,
            dplyr::across(dplyr::all_of(cols), ~ list(cell[.x])),
            conf_quantile = mean(conf_quantile)
        )
    } else {
        dplyr::summarise(
            res,
            dplyr::across(dplyr::all_of(cols), ~ list(cell[.x]))
        )
    }
    dplyr::mutate(out, cp_alpha = alpha)
}

#' Check argument's validity passed to conformeR
#' @noRd
check_arg <- function(ok, msg) if (!ok) stop(msg, call. = FALSE)
