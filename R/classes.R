#' Guide-level log-fold-change container
#'
#' @description
#' Internal S4 container used to store guide-level log-fold-change arrays and
#' their dimension metadata.
#'
#' @slot data Numeric array of guide-level values.
#' @slot space Character vector naming biological dimensions.
#' @slot replicates Character vector naming replicate dimensions.
#'
#' @keywords internal
#' @exportClass gRNA_LFC

setClass(
  "gRNA_LFC",
  slots = list(
    "data" = "array",
    "space" = "character",
    "replicates" = "character"
  )
)

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
#' @keywords internal
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

#' Base class for CRISPR genetic-interaction screens
#'
#' @description
#' `ScreenBase` stores all intermediate and final data used by the CeRberus GI
#' aggregation workflow. Concrete screen types inherit from this class.
#'
#' @slot guideLFCs A [`gRNA_LFC-class`] object.
#' @slot guideGIs A [`gRNA_GI-class`] object.
#' @slot limma_models List of fitted limma model objects.
#' @slot geneGIs Numeric array of aggregated gene-level GI scores, p-values, and
#'   FDR values.
#' @slot screen_attr List of inferred screen attributes.
#' @slot dupCorrelation Numeric duplicate-correlation estimate(s).
#' @slot metadata List of run metadata and standardized input data.
#' @slot checks List of screen-design checks.
#' @slot errors List of recoverable modelling errors.
#'
#' @keywords internal
#' @exportClass ScreenBase

setClass(
  "ScreenBase",
  slots = list(
    "guideLFCs" = "gRNA_LFC",
    "guideGIs" = "gRNA_GI",
    "limma_models" = "list",
    "geneGIs" = "array",
    "screen_attr" = "list",
    #"blocks" = "character",
    #"block_description" = "list",
    "dupCorrelation" = "numeric",
    "metadata" = "list",
    "checks" = "list",
    "errors" = "list"
  )
)

#' Fixed-pair CRISPR screen
#'
#' @description
#' S4 class for screens where each guide-pair maps to a fixed gene pair.
#'
#' @keywords internal
#' @exportClass FixedPairScreen

setClass("FixedPairScreen", contains = "ScreenBase")

#' Multiplex CRISPR screen
#'
#' @description
#' S4 class for multiplex screens where query genes are modelled against a
#' library-gene axis.
#'
#' @keywords internal
#' @exportClass MultiplexScreen

setClass("MultiplexScreen", contains = "ScreenBase")

#' Position-agnostic multiplex CRISPR screen
#'
#' @description
#' S4 class for multiplex screens whose directional guide-level scores have
#' been symmetrized before gene-pair aggregation.
#'
#' @slot symmGeneGIs Data table with one row per unordered gene pair.
#'
#' @keywords internal
#' @exportClass PosAgnMultiplexScreen

setClass(
  "PosAgnMultiplexScreen",
  contains = "MultiplexScreen",
  slots = list("symmGeneGIs" = "data.table")
)


# Need to be called without a loop to work

setAs(from = "ScreenBase", to = "FixedPairScreen", function(from) {
  obj <- new("FixedPairScreen")
  for (s in slotNames("ScreenBase")) {
    slot(obj, s) <- slot(from, s)
  }
  class(obj) <- "FixedPairScreen"
  return(obj)
})

setAs(from = "ScreenBase", to = "MultiplexScreen", function(from) {
  obj <- new("MultiplexScreen")
  for (s in slotNames("ScreenBase")) {
    slot(obj, s) <- slot(from, s)
  }
  class(obj) <- "MultiplexScreen"
  return(obj)
})

setAs(from = "ScreenBase", to = "PosAgnMultiplexScreen", function(from) {
  obj <- new("PosAgnMultiplexScreen")
  for (s in slotNames("ScreenBase")) {
    slot(obj, s) <- slot(from, s)
  }
  class(obj) <- "PosAgnMultiplexScreen"
  return(obj)
})
