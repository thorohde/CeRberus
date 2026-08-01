#' Guide-level genetic-interaction container
#'
#' @description
#' Internal S4 container used by CeRberus to store guide-level genetic
#' interaction scores, replicate metadata, and duplicate-correlation blocking
#' information.
#'
#' @slot data Numeric array of guide-level GI scores.
#' @slot space Character vector naming biological dimensions, for example
#'   `"gene_pair"` or `c("query_gene", "library_gene")`.
#' @slot replicates Character vector naming available replicate dimensions.
#' @slot block_layer Character scalar naming the replicate layer used for limma
#'   blocking, or `character(0)` if none is used.
#' @slot blocks Character vector of block assignments passed to limma.
#' @slot use_blocks Logical scalar indicating whether blocking is active.
#' @slot block_description Character vector describing flattened replicate
#'   columns.
#' @slot collapse Character vector of replicate layers collapsed before model
#'   fitting.
#'
#' @exportClass gRNA_GI

setClass(
  "gRNA_GI",
  slots = list(
    "data" = "array",
    "space" = "character",
    "replicates" = "character",
    "block_layer" = "character",
    "blocks" = "character",
    "use_blocks" = "logical",
    "block_description" = "character",
    "collapse" = "character"
  )
)

setValidity("gRNA_GI", function(object) {
  structure_validity <- .validate_guide_container_structure(object)
  if (!isTRUE(structure_validity)) {
    return(structure_validity)
  }

  if (
    length(object@block_layer) > 1L ||
      anyNA(object@block_layer) ||
      any(!nzchar(object@block_layer))
  ) {
    return("'block_layer' must be empty or one non-empty, non-missing name.")
  }

  if (any(object@block_layer %in% object@space)) {
    return("A biological 'space' dimension cannot be used as 'block_layer'.")
  }

  if (length(object@use_blocks) != 1L || is.na(object@use_blocks)) {
    return("'use_blocks' must be TRUE or FALSE.")
  }

  if (isTRUE(object@use_blocks) && length(object@blocks) == 0L) {
    return("'blocks' must be populated when 'use_blocks' is TRUE.")
  }

  if (
    length(object@block_description) > 0L &&
      length(object@data) > 0L &&
      length(object@block_description) !=
        dim(object@data)[[length(dim(object@data))]]
  ) {
    return("'block_description' must describe the final 'data' dimension.")
  }

  if (
    isTRUE(object@use_blocks) &&
      length(object@block_description) > 0L &&
      length(object@blocks) != length(object@block_description)
  ) {
    return("'blocks' and 'block_description' must have equal lengths.")
  }

  TRUE
})
