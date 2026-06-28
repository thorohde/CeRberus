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
