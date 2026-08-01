fill_gRNA_LFCs <- function(gRNA_LFC, input, value_var = "LFC") {
  fill_gRNA_GIs(gRNA_LFC, input, value_var = value_var)
}

#####

#' Compute positional gene main effects from guide-level LFCs
#'
#' @description
#' Computes one main effect per query gene and one main effect per library gene
#' from a directional guide-level log-fold-change array. Query main effects are
#' means over all library genes and retained observations; library main effects
#' are means over all query genes and retained observations.
#'
#' Missing values are removed when calculating means. A gene with only missing
#' measurements receives `NA_real_`. The guide-level LFC array is not modified.
#'
#' @param gRNA_LFC A [`gRNA_LFC-class`] object containing a directional LFC
#'   array. The first two dimensions must be named `query_gene` and
#'   `library_gene`.
#'
#' @return The input `gRNA_LFC` object with named `query_main_effects` and
#'   `library_main_effects` slots populated.
#'
#' @export
compute_gene_main_effects <- function(gRNA_LFC) {
  stopifnot(
    "gRNA_LFC must be a gRNA_LFC object." = methods::is(gRNA_LFC, "gRNA_LFC")
  )
  stopifnot(
    "The LFC array must contain at least two dimensions." = length(
      dim(gRNA_LFC@data)
    ) >=
      2L
  )

  query_genes <- dimnames(gRNA_LFC@data)[[1L]]
  library_genes <- dimnames(gRNA_LFC@data)[[2L]]

  stopifnot(
    "The first two LFC dimensions must be query_gene and library_gene." = identical(
      names(dimnames(gRNA_LFC@data))[1:2],
      c("query_gene", "library_gene")
    ),
    "The query-gene dimension must have non-empty names." = !is.null(
      query_genes
    ) &&
      !anyNA(query_genes) &&
      all(nzchar(query_genes)),
    "The library-gene dimension must have non-empty names." = !is.null(
      library_genes
    ) &&
      !anyNA(library_genes) &&
      all(nzchar(library_genes))
  )

  query_effects <- apply(
    X = gRNA_LFC@data,
    MARGIN = 1L,
    FUN = mean,
    na.rm = TRUE
  )
  library_effects <- apply(
    X = gRNA_LFC@data,
    MARGIN = 2L,
    FUN = mean,
    na.rm = TRUE
  )

  query_effects[is.nan(query_effects)] <- NA_real_
  library_effects[is.nan(library_effects)] <- NA_real_

  gRNA_LFC@query_main_effects <- stats::setNames(
    as.numeric(query_effects),
    query_genes
  )
  gRNA_LFC@library_main_effects <- stats::setNames(
    as.numeric(library_effects),
    library_genes
  )

  methods::validObject(gRNA_LFC)
  gRNA_LFC
}

#####

#' @describeIn collapse_replicates Collapse replicate dimensions in a
#'   [`gRNA_LFC-class`] object by averaging over layers listed in the `collapse`
#'   slot.
#' @aliases collapse_replicates,gRNA_LFC-method
setMethod(
  "collapse_replicates",
  signature = signature(.x = "gRNA_LFC"),
  function(.x) {
    if (length(.x@collapse) == 0L) {
      return(.x)
    }

    current_dims <- c(.x@space, .x@replicates)

    if (!all(.x@collapse %in% .x@replicates)) {
      stop(
        "All collapse layers must be replicate dimensions. Invalid layer(s): ",
        paste(setdiff(.x@collapse, .x@replicates), collapse = ", "),
        call. = FALSE
      )
    }

    keep_dims <- setdiff(current_dims, .x@collapse)
    keep_margin <- match(keep_dims, current_dims)
    keep_dimnames <- dimnames(.x@data)[keep_margin]
    names(keep_dimnames) <- keep_dims

    .x@data <- array(
      data = apply(
        X = .x@data,
        MARGIN = keep_margin,
        FUN = mean,
        na.rm = TRUE
      ),
      dim = purrr::map_int(keep_dimnames, length),
      dimnames = keep_dimnames
    )
    .x@replicates <- setdiff(.x@replicates, .x@collapse)

    return(.x)
  }
)

#####

flatten_guide_lfcs <- function(.x) {
  flatten_formula <- stats::as.formula(paste0(
    paste0(.x@space, collapse = " ~ "),
    " ~ ",
    paste0(.x@replicates, collapse = " + ")
  ))

  .x@data <- .x@data |>
    flatten_array(c(.x@space, .x@replicates)) |>
    reshape2::acast(flatten_formula)

  return(.x)
}
