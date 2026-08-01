.validate_guide_container_structure <- function(object) {
  dimension_fields <- c(object@space, object@replicates)

  if (anyNA(dimension_fields) || any(!nzchar(dimension_fields))) {
    return(
      "'space' and 'replicates' must contain non-empty, non-missing names."
    )
  }

  if (anyDuplicated(dimension_fields)) {
    return("'space' and 'replicates' must contain unique dimension names.")
  }

  if (
    anyNA(object@collapse) ||
      any(!nzchar(object@collapse)) ||
      anyDuplicated(object@collapse)
  ) {
    return(
      "'collapse' must contain unique, non-empty, non-missing layer names."
    )
  }

  if (any(object@collapse %in% object@space)) {
    return("Biological 'space' dimensions cannot be collapsed as replicates.")
  }

  if (length(object@data) == 0L) {
    return(TRUE)
  }

  data_dims <- dim(object@data)
  expected_ranks <- unique(c(
    length(dimension_fields),
    length(object@space) + as.integer(length(object@replicates) > 0L)
  ))

  if (!length(data_dims) %in% expected_ranks) {
    return(paste0(
      "The 'data' array rank must match the biological and replicate metadata ",
      "before or after replicate flattening."
    ))
  }

  data_dimnames <- dimnames(object@data)
  if (!is.null(data_dimnames)) {
    named_dimensions <- names(data_dimnames)

    if (
      !is.null(named_dimensions) &&
        any(nzchar(named_dimensions)) &&
        !identical(named_dimensions[seq_along(object@space)], object@space)
    ) {
      return("Leading 'data' dimension names must match the 'space' metadata.")
    }

    if (any(purrr::map_lgl(data_dimnames, anyDuplicated))) {
      return("Each 'data' dimension must have unique labels.")
    }
  }

  TRUE
}

#' Guide-level log-fold-change container
#'
#' @description
#' Internal S4 container used to store guide-level log-fold-change arrays and
#' their dimension metadata.
#'
#' @slot data Numeric array of guide-level values.
#' @slot space Character vector naming biological dimensions.
#' @slot replicates Character vector naming replicate dimensions.
#' @slot collapse Character vector of replicate layers collapsed during
#'   construction.
#' @slot query_main_effects Named numeric vector of query-gene main effects.
#' @slot library_main_effects Named numeric vector of library-gene main effects.
#'
#' @exportClass gRNA_LFC

setClass(
  "gRNA_LFC",
  slots = list(
    "data" = "array",
    "space" = "character",
    "replicates" = "character",
    "collapse" = "character",
    "query_main_effects" = "numeric",
    "library_main_effects" = "numeric"
  )
)

setValidity("gRNA_LFC", function(object) {
  structure_validity <- .validate_guide_container_structure(object)
  if (!isTRUE(structure_validity)) {
    return(structure_validity)
  }

  main_effect_slots <- c("query_main_effects", "library_main_effects")

  for (.slot in main_effect_slots) {
    .values <- methods::slot(object, .slot)
    .names <- names(.values)

    if (
      length(.values) > 0L &&
        (is.null(.names) ||
          anyNA(.names) ||
          any(!nzchar(.names)) ||
          anyDuplicated(.names))
    ) {
      return(paste0(
        "'",
        .slot,
        "' must be a numeric vector with unique, non-empty names."
      ))
    }
  }

  TRUE
})
