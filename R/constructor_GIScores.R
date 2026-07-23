#' Construct a genetic-interaction screen object
#'
#' @description
#' `GIScores()` imports guide-level genetic-interaction scores, infers the
#' screen structure, stores the scores in CeRberus S4 containers, and prepares
#' the object for duplicate-correlation estimation and limma-based aggregation.
#'
#' The input table is expected in long format with one guide-pair observation
#' per row. Column names can be customized with the `*_col` arguments and are
#' standardized internally.
#'
#' @param input A data frame or data.table containing guide-level GI scores.
#' @param query_col Name of the column containing query-gene identifiers.
#' @param lib_col Name of the column containing library-gene identifiers.
#' @param bio_rep_col Name of the biological-replicate column.
#' @param tech_rep_col Name of the technical-replicate column.
#' @param guide_col Name of the guide-pair column.
#' @param gi_col Name of the numeric genetic-interaction score column.
#' @param collapse_layers Optional character vector of replicate layers to
#'   average before model fitting, for example `"guide_pair"`, `"tech_rep"`, or
#'   `"bio_rep"`.
#' @param block_layer Optional replicate layer used as the limma blocking layer
#'   for duplicate-correlation modelling.
#' @param force_fixed_pair Deprecated logical override. Use
#'   `screen_type = "fixed_pair"` instead.
#' @param screen_type Character scalar controlling screen classification.
#'   `"auto"` uses the inferred screen attributes, while `"fixed_pair"` and
#'   `"multiplex"` explicitly select the corresponding screen structure.
#' @param pos_agnostic Logical. If `TRUE`, opt into position-agnostic analysis
#'   for an inferred multiplex screen by averaging both orientations of each
#'   gene pair. The default, `FALSE`, retains directional query-by-library
#'   analysis.
#' @param symmetric_analysis_method Character scalar specifying how
#'   position-agnostic multiplex screens are modeled. `preaverage` averages
#'   both pair orientations before fitting one limma model per query gene.
#'   `global_preaverage` performs the same orientation averaging, converts the
#'   result to one row per unordered gene pair, and fits one limma model across
#'   all unordered pairs. This argument is only used when
#'   `pos_agnostic = TRUE`; directional analysis remains the default.
#' @param verbose Logical. If `TRUE`, print a short screen summary.
#'
#' @return An S4 object inheriting from `ScreenBase`, typically a
#'   `FixedPairScreen`, `MultiplexScreen`, or `PosAgnMultiplexScreen`.
#'
#' @export

GIScores <- function(
  input,
  query_col = "query_gene",
  lib_col = "library_gene",
  bio_rep_col = "bio_rep",
  tech_rep_col = "tech_rep",
  guide_col = "guide_pair",
  gi_col = "GI",
  collapse_layers = NULL,
  block_layer = NULL,
  force_fixed_pair = FALSE,
  pos_agnostic = FALSE,
  symmetric_analysis_method = "preaverage",
  verbose = FALSE,
  screen_type = c("auto", "fixed_pair", "multiplex")
) {
  screen_type <- match.arg(screen_type)

  if (
    !is.logical(force_fixed_pair) ||
      length(force_fixed_pair) != 1L ||
      is.na(force_fixed_pair)
  ) {
    stop("force_fixed_pair must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(force_fixed_pair)) {
    .Deprecated(
      msg = paste0(
        "force_fixed_pair is deprecated; use ",
        "screen_type = \"fixed_pair\" instead."
      )
    )

    if (!identical(screen_type, "auto")) {
      stop(
        "Do not use force_fixed_pair together with a non-'auto' screen_type.",
        call. = FALSE
      )
    }

    screen_type <- "fixed_pair"
  }

  if (
    !is.logical(pos_agnostic) ||
      length(pos_agnostic) != 1L ||
      is.na(pos_agnostic)
  ) {
    stop(
      "pos_agnostic must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (isTRUE(pos_agnostic)) {
    symmetric_analysis_method <- validate_symmetric_analysis_method(
      symmetric_analysis_method
    )
  }

  gi_obj <- new(
    "ScreenBase",
    metadata = list(
      input = input,
      query_col = query_col,
      lib_col = lib_col,
      bio_rep_col = bio_rep_col,
      tech_rep_col = tech_rep_col,
      guide_col = guide_col,
      gi_col = gi_col,
      requested_screen_type = screen_type,
      symmetric_analysis_method = symmetric_analysis_method
    )
  )

  gi_obj@guideGIs@collapse <- as.character(collapse_layers)
  gi_obj@guideGIs@block_layer <- as.character(block_layer)

  gi_obj <- import_scores(gi_obj)
  gi_obj <- get_screen_attributes(gi_obj)
  gi_obj <- run_checks(gi_obj)
  gi_obj <- set_screen_type(gi_obj)

  if (pos_agnostic && !methods::is(gi_obj, "MultiplexScreen")) {
    stop(
      "Position-agnostic analysis is only available for inferred multiplex screens.",
      call. = FALSE
    )
  }

  gi_obj@guideGIs <- gi_obj@guideGIs |> fill_gRNA_GIs(gi_obj@metadata$input)
  gi_obj@guideGIs <- gi_obj@guideGIs |> collapse_replicates()
  gi_obj@guideGIs <- gi_obj@guideGIs |> flatten_guide_gis()

  if (methods::is(gi_obj, "MultiplexScreen") && pos_agnostic) {
    query_genes <- rownames(gi_obj@guideGIs@data)
    library_genes <- colnames(gi_obj@guideGIs@data)

    if (!setequal(query_genes, library_genes)) {
      stop(
        "Position-agnostic analysis requires identical query and library gene sets.",
        call. = FALSE
      )
    }

    gi_obj@guideGIs@data <- gi_obj@guideGIs@data[
      query_genes,
      query_genes,
      ,
      drop = FALSE
    ]

    for (.r in gi_obj@guideGIs@block_description) {
      gi_obj@guideGIs@data[,, .r] <- make_symmetric(
        gi_obj@guideGIs@data[,, .r]
      )
    }

    if (identical(symmetric_analysis_method, "global_preaverage")) {
      gi_obj@guideGIs@data <- flatten_symmetric_pairs(
        .arr = gi_obj@guideGIs@data,
        pairs = gi_obj@screen_attr$unique_pairs
      )

      gi_obj@guideGIs@space <- "gene_pair"
    }

    gi_obj <- methods::as(
      object = gi_obj,
      Class = "PosAgnMultiplexScreen"
    )
  }

  if (isTRUE(verbose)) {
    screen_report(gi_obj)
  }

  return(gi_obj)
}
