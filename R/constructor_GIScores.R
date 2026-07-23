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
#' @param force_fixed_pair Logical. If `TRUE`, force construction as a
#'   fixed-pair screen regardless of the inferred screen attributes.
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
  verbose = FALSE
) {
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

  GI_obj <- new(
    "ScreenBase",
    metadata = list(
      input = input,
      query_col = query_col,
      lib_col = lib_col,
      bio_rep_col = bio_rep_col,
      tech_rep_col = tech_rep_col,
      guide_col = guide_col,
      gi_col = gi_col,
      force_fixed_pair = force_fixed_pair,
      symmetric_analysis_method = symmetric_analysis_method
    )
  )

  GI_obj@guideGIs@collapse <- as.character(collapse_layers)
  GI_obj@guideGIs@block_layer <- as.character(block_layer)

  GI_obj <- import_scores(GI_obj)
  GI_obj <- get_screen_attributes(GI_obj)
  GI_obj <- run_checks(GI_obj)
  GI_obj <- set_screenType(GI_obj)

  if (pos_agnostic && !methods::is(GI_obj, "MultiplexScreen")) {
    stop(
      "Position-agnostic analysis is only available for inferred multiplex screens.",
      call. = FALSE
    )
  }

  GI_obj@guideGIs <- GI_obj@guideGIs |> fill_gRNA_GIs(GI_obj@metadata$input)
  GI_obj@guideGIs <- GI_obj@guideGIs |> collapse_replicates()
  GI_obj@guideGIs <- GI_obj@guideGIs |> flatten_guideGIs()

  if (methods::is(GI_obj, "MultiplexScreen") && pos_agnostic) {
    query_genes <- rownames(GI_obj@guideGIs@data)
    library_genes <- colnames(GI_obj@guideGIs@data)

    if (!setequal(query_genes, library_genes)) {
      stop(
        "Position-agnostic analysis requires identical query and library gene sets.",
        call. = FALSE
      )
    }

    GI_obj@guideGIs@data <- GI_obj@guideGIs@data[
      query_genes,
      query_genes,
      ,
      drop = FALSE
    ]

    for (.r in GI_obj@guideGIs@block_description) {
      GI_obj@guideGIs@data[,, .r] <- makeSymmetric(
        GI_obj@guideGIs@data[,, .r]
      )
    }

    if (identical(symmetric_analysis_method, "global_preaverage")) {
      GI_obj@guideGIs@data <- flatten_symmetric_pairs(
        .arr = GI_obj@guideGIs@data,
        pairs = GI_obj@screen_attr$unique_pairs
      )

      GI_obj@guideGIs@space <- "gene_pair"
    }

    GI_obj <- methods::as(
      object = GI_obj,
      Class = "PosAgnMultiplexScreen"
    )
  }

  if (isTRUE(verbose)) {
    screenReport(GI_obj)
  }

  ##### rework?
  #GI_obj <- compute_dupCorrelation(GI_obj)

  #GI_obj@metadata$used_block <- block_layer

  return(GI_obj)
}
