#####

symmetric_analysis_methods <- function() {
  c("preaverage", "global_preaverage")
}

validate_symmetric_analysis_method <- function(
  symmetric_analysis_method = "preaverage"
) {
  if (
    !is.character(symmetric_analysis_method) ||
      length(symmetric_analysis_method) != 1L ||
      is.na(symmetric_analysis_method)
  ) {
    stop(
      "symmetric_analysis_method must be a single non-missing character value.",
      call. = FALSE
    )
  }

  match.arg(
    symmetric_analysis_method,
    symmetric_analysis_methods()
  )
}

#####

flatten_symmetric_pairs <- function(.arr, pairs, sep = ";") {
  stopifnot(
    "The symmetric guide-level input must be a three-dimensional array." = length(
      dim(.arr)
    ) ==
      3L,
    "pairs must be a character vector." = is.character(pairs),
    "pairs must not contain missing values." = !anyNA(pairs),
    "pairs must be unique." = !anyDuplicated(pairs)
  )

  genes1 <- stringr::str_split_i(pairs, sep, 1)
  genes2 <- stringr::str_split_i(pairs, sep, 2)

  stopifnot(
    "Some first-position genes are absent from the query-gene dimension." = all(
      genes1 %in% rownames(.arr)
    ),
    "Some second-position genes are absent from the library-gene dimension." = all(
      genes2 %in% colnames(.arr)
    )
  )

  n_observations <- dim(.arr)[3L]

  output <- base::matrix(
    data = base::unlist(
      lapply(seq_along(pairs), function(.i) {
        as.numeric(.arr[genes1[.i], genes2[.i], , drop = TRUE])
      }),
      use.names = FALSE
    ),
    nrow = length(pairs),
    ncol = n_observations,
    byrow = TRUE,
    dimnames = list(
      gene_pair = pairs,
      replicate = dimnames(.arr)[[3L]]
    )
  )

  # mean(..., na.rm = TRUE) returns NaN when both orientations are missing.
  # limma should receive those values as ordinary missing observations.
  output[is.nan(output)] <- NA_real_

  return(output)
}

#####

get_symmetric_analysis_method <- function(GI_obj) {
  method <- GI_obj@metadata$symmetric_analysis_method

  if (is.null(method) || length(method) == 0L) {
    return("preaverage")
  }

  validate_symmetric_analysis_method(method)
}

#####

makeSymmetric <- function(.x) {
  output <- base::apply(
    abind::abind(.x, base::t(.x), along = 3),
    1:2,
    base::mean,
    na.rm = TRUE
  )

  output[is.nan(output)] <- NA_real_

  return(output)
}

#####

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

#####

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

#####

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

#####
