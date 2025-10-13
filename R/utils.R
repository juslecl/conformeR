# Subset SCE object according to a condition
subset_by_condition <- function(sce, condition, value) {
  sce[, colData(sce)[[condition]] == value]
}


# helper: convert SCE to wide tibble with conf_group
sce_to_wide_tibble <- function(sce,obs_condition) {
  meta_df <- as.data.frame(colData(sce)) %>%
    select(conf_group, all_of(obs_condition))

  expr_df <- as.data.frame(t(assay(sce)))  # transpose: rows=cells, cols=genes
  colnames(expr_df) <- rownames(sce)

  bind_cols(meta_df, expr_df) %>% as_tibble()
}

# Derives p-values from confidence intervals
interval_to_pval <- function(intervals_df){
  pval_df <-intervals_df |>
    group_by(cell_id,conf_group,gene) |>
    summarise(
      pvalue = (1 + sum(covered)) / (1 + n()),
      .groups = "drop"
    )
  return(pval_df)
}

# Maps p-values to fdr
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
  fdr_val <- ifelse(sum(pval < cutoff) == 0, 1, mean(lfdr[pval < cutoff]))
  Rg <- sum(pval <= cutoff)
  int <- int |> mutate(Rg=Rg,fdr=fdr_val)
  return(int)
}

# computes combined FDR
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
        mean(Rg * fdr, na.rm = TRUE),
        sum(Rg * fdr, na.rm = TRUE) / sum(Rg, na.rm = TRUE)
      ),
      comb_fdr = if_else(Rg == 0, 1, comb_fdr)
    )
  fdr_tab
}
