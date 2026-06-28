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
