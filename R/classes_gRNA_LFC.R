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
#' @exportClass gRNA_LFC

setClass(
  "gRNA_LFC",
  slots = list(
    "data" = "array",
    "space" = "character",
    "replicates" = "character"
  )
)
