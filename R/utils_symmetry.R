makeSymmetric <- function(.x) {
  base::apply(
    abind::abind(.x, base::t(.x), along = 3),
    1:2,
    base::mean,
    na.rm = TRUE
  )
}


gather_symmetric_scores <- function(pairs, .arr, sep = ";") {
  #  if (!isSymmetric(.arr)) {
  #    warning(stringr::str_c("Input array is asymmetric! ABBA: ",
  #                           round(abba_cor(.arr), 3)))}

  genes1 <- str_split_i(pairs, sep, 1)
  genes2 <- str_split_i(pairs, sep, 2)

  stopifnot(
    length(setdiff(genes1, rownames(.arr))) == 0,
    length(setdiff(genes2, colnames(.arr))) == 0
  )

  return(purrr::map2_dbl(genes1, genes2, \(.g1, .g2) {
    .arr[.g1, .g2]
  }))
}


#' Compute balanced FDR values for position-agnostic gene pairs
#'
#' @description
#' For each unordered gene pair, `balanced_FDR()` builds a local set of
#' directional p-values involving the two genes and adjusts those p-values for
#' multiple testing. The returned value for each pair is the adjusted p-value
#' corresponding to the original pair orientation.
#'
#' This helper is used when aggregating symmetric or position-agnostic multiplex
#' screens, where both `A;B` and `B;A` directional tests may contribute to the
#' final gene-pair-level result.
#'
#' @param pairs Character vector of gene-pair identifiers in the form
#'   `"gene1;gene2"`.
#' @param pval_array Numeric matrix-like object of p-values with query genes in
#'   rows and library genes in columns.
#' @param fdr_method Multiple-testing correction method passed to
#'   [stats::p.adjust()].
#'
#' @return A numeric vector of adjusted p-values with one value per element of
#'   `pairs`.
#'

balanced_FDR <- function(pairs, pval_array, fdr_method) {
  pair_template <- \(pair) {
    .g1 <- str_split_i(pair, ";", 1)
    .g2 <- str_split_i(pair, ";", 2)
    .d <- data.table(
      gene_pair = pair,
      g1 = c(rep(.g1, ncol(pval_array)), setdiff(rownames(pval_array), .g1)),
      g2 = c(colnames(pval_array), rep(.g2, nrow(pval_array) - 1))
    )
    .d[, gene_pair2 := stringr::str_c(g1, ";", g2)]
    .d
  }

  return(purrr::map_dbl(pairs, \(pair) {
    .d <- pair_template(pair)
    .d[, pval := purrr::map2_dbl(g1, g2, ~ pval_array[.x, .y])]
    .d[, FDR := stats::p.adjust(pval, method = fdr_method)]
    .d[gene_pair == gene_pair2, get("FDR")]
  }))
}
