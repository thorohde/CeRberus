#' Internal screen-design container
#'
#' @description
#' `ScreenDesign` stores the canonical gene sets and pair structure inferred
#' from a genetic-interaction screen. It is an implementation detail of
#' [`ScreenBase-class`] objects.
#'
#' @slot contrasts Optional contrast specification used by downstream models.
#' @slot query_genes Unique query-gene names in canonical input order.
#' @slot library_genes Unique library-gene names in canonical input order.
#' @slot all_genes Ordered union of query and library genes.
#' @slot query_genes_not_in_lib Query genes absent from the library axis.
#' @slot library_genes_not_in_query Library genes absent from the query axis.
#' @slot n_query_genes Number of query genes.
#' @slot n_lib_genes Number of library genes.
#' @slot n_all_genes Number of genes across both axes.
#' @slot observations_per_query Named non-negative integer vector with one
#'   observation count per query gene in canonical order.
#' @slot all_pairs Unique directional gene pairs encoded as `"query;library"`.
#' @slot unique_pairs Canonical unordered forms of `all_pairs` in first-seen
#'   order.
#'
#' @keywords internal

setClass(
  "ScreenDesign",
  slots = list(
    "contrasts" = "ANY",
    "query_genes" = "character",
    "library_genes" = "character",
    "all_genes" = "character",
    "query_genes_not_in_lib" = "character",
    "library_genes_not_in_query" = "character",
    "n_query_genes" = "integer",
    "n_lib_genes" = "integer",
    "n_all_genes" = "integer",
    "observations_per_query" = "integer",
    "all_pairs" = "character",
    "unique_pairs" = "character"
  ),
  prototype = list(
    contrasts = NULL,
    query_genes = character(),
    library_genes = character(),
    all_genes = character(),
    query_genes_not_in_lib = character(),
    library_genes_not_in_query = character(),
    n_query_genes = 0L,
    n_lib_genes = 0L,
    n_all_genes = 0L,
    observations_per_query = integer(),
    all_pairs = character(),
    unique_pairs = character()
  )
)

setValidity("ScreenDesign", function(object) {
  gene_fields <- c(
    "query_genes",
    "library_genes",
    "all_genes",
    "query_genes_not_in_lib",
    "library_genes_not_in_query"
  )

  for (.field in gene_fields) {
    .genes <- methods::slot(object, .field)
    if (anyNA(.genes) || any(!nzchar(.genes)) || anyDuplicated(.genes)) {
      return(paste0(
        "'",
        .field,
        "' must contain unique, non-empty, non-missing gene names."
      ))
    }
  }

  count_map <- c(
    n_query_genes = "query_genes",
    n_lib_genes = "library_genes",
    n_all_genes = "all_genes"
  )

  for (.count_field in names(count_map)) {
    .count <- methods::slot(object, .count_field)
    .gene_field <- count_map[[.count_field]]

    if (
      length(.count) != 1L ||
        is.na(.count) ||
        .count < 0L ||
        .count != length(methods::slot(object, .gene_field))
    ) {
      return(paste0(
        "'",
        .count_field,
        "' must equal the length of '",
        .gene_field,
        "'."
      ))
    }
  }

  expected_all_genes <- union(object@query_genes, object@library_genes)
  if (!identical(object@all_genes, expected_all_genes)) {
    return("'all_genes' must equal union(query_genes, library_genes).")
  }

  expected_query_only <- setdiff(object@query_genes, object@library_genes)
  if (!identical(object@query_genes_not_in_lib, expected_query_only)) {
    return(
      "'query_genes_not_in_lib' must equal setdiff(query_genes, library_genes)."
    )
  }

  expected_library_only <- setdiff(object@library_genes, object@query_genes)
  if (!identical(object@library_genes_not_in_query, expected_library_only)) {
    return(paste0(
      "'library_genes_not_in_query' must equal ",
      "setdiff(library_genes, query_genes)."
    ))
  }

  observations <- object@observations_per_query
  if (
    length(observations) != length(object@query_genes) ||
      !identical(names(observations), object@query_genes) ||
      anyNA(observations) ||
      any(observations < 0L)
  ) {
    return(paste0(
      "'observations_per_query' must be a named, non-negative integer vector ",
      "with one value per query gene in canonical order."
    ))
  }

  pair_fields <- c("all_pairs", "unique_pairs")
  for (.field in pair_fields) {
    .pairs <- methods::slot(object, .field)
    if (anyNA(.pairs) || any(!nzchar(.pairs)) || anyDuplicated(.pairs)) {
      return(paste0(
        "'",
        .field,
        "' must contain unique, non-empty, non-missing gene pairs."
      ))
    }
  }

  if (length(object@all_pairs) > 0L) {
    pair_parts <- strsplit(object@all_pairs, ";", fixed = TRUE)
    if (
      any(lengths(pair_parts) != 2L) ||
        any(!unlist(pair_parts) %in% object@all_genes)
    ) {
      return("'all_pairs' must contain ';'-separated genes from 'all_genes'.")
    }

    expected_unique_pairs <- unique(sort_gene_pairs(
      purrr::map_chr(pair_parts, 1L),
      purrr::map_chr(pair_parts, 2L)
    ))
    if (!identical(object@unique_pairs, expected_unique_pairs)) {
      return(
        "'unique_pairs' must be the canonical unordered form of 'all_pairs'."
      )
    }
  } else if (length(object@unique_pairs) > 0L) {
    return("'unique_pairs' must be empty when 'all_pairs' is empty.")
  }

  TRUE
})
