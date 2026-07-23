#' Base class for CRISPR genetic-interaction screens
#'
#' @description
#' `ScreenBase` stores all intermediate and final data used by the CeRberus GI
#' aggregation workflow. Concrete screen types inherit from this class.
#'
#' @slot guideLFCs A [`gRNA_LFC-class`] object.
#' @slot guideGIs A [`gRNA_GI-class`] object.
#' @slot limma_models Fitted limma output. Multiplex analyses usually store a
#'   list of per-query models; global position-agnostic analysis stores one
#'   `MArrayLM` object covering all unordered pairs.
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
#' S4 class for multiplex screens analyzed without distinguishing query-library
#' orientation in the final gene-pair output. Pair orientations are averaged
#' before model fitting. Depending on the selected symmetric analysis method,
#' limma models are fitted either per query gene or once across all unordered
#' gene pairs.
#'
#' @slot symmGeneGIs Data table with one row per unordered gene pair.
#'
#' @keywords internal
#' @exportClass PosAgnMultiplexScreen

setClass(
  "PosAgnMultiplexScreen",
  contains = "MultiplexScreen",
  slots = list("symmGeneGIs" = "data.table"),
  prototype = list(symmGeneGIs = data.table::data.table())
)


# Need to be called without a loop to work

setAs(from = "ScreenBase", to = "FixedPairScreen", function(from) {
  obj <- new("FixedPairScreen")
  for (s in slotNames("ScreenBase")) {
    slot(obj, s) <- slot(from, s)
  }
  return(obj)
})

setAs(from = "ScreenBase", to = "MultiplexScreen", function(from) {
  obj <- new("MultiplexScreen")
  for (s in slotNames("ScreenBase")) {
    slot(obj, s) <- slot(from, s)
  }
  return(obj)
})

setAs(from = "ScreenBase", to = "PosAgnMultiplexScreen", function(from) {
  obj <- new("PosAgnMultiplexScreen")
  for (s in slotNames("ScreenBase")) {
    slot(obj, s) <- slot(from, s)
  }
  return(obj)
})
