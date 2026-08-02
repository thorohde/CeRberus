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
#' For screens with guide LFC input, `GIScores()` automatically computes
#' separate query- and library-position gene main effects. These are retained as
#' provenance in the `gRNA_LFC` container. After model fitting and
#' `collect_gis()`, `gi_df()` adds `query_main_effect` and
#' `library_main_effect` columns to directional results. Position-agnostic
#' results instead contain `gene1_main_effect` and `gene2_main_effect`, based on
#' the equal mean of each gene's query and library effects. These columns can be
#' used to construct covariates for downstream procedures such as IHW; CeRberus
#' continues to apply the requested standard `p.adjust()` method internally.
#'
#' @param input A data frame or data.table containing guide-level GI scores.
#' @param query_col Name of the column containing query-gene identifiers.
#' @param lib_col Name of the column containing library-gene identifiers.
#' @param bio_rep_col Name of the biological-replicate column.
#' @param tech_rep_col Name of the technical-replicate column.
#' @param guide_col Name of the guide-pair column.
#' @param gi_col Name of the numeric genetic-interaction score column.
#' @param lfc_col Optional name of the numeric guide log-fold-change column.
#'   If `NULL`, a column named `LFC` is used when present. If no such column is
#'   present, the returned object's guide-level LFC container remains empty.
#'   Supplied LFCs also trigger automatic positional gene main-effect
#'   computation for every screen design.
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
  screen_type = c("auto", "fixed_pair", "multiplex"),
  lfc_col = NULL
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
      lfc_col = lfc_col,
      requested_screen_type = screen_type,
      symmetric_analysis_method = symmetric_analysis_method
    )
  )

  gi_obj@guideGIs@collapse <- as.character(collapse_layers)
  gi_obj@guideGIs@block_layer <- as.character(block_layer)
  gi_obj@guideLFCs@collapse <- as.character(collapse_layers)

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

  if ("LFC" %in% colnames(gi_obj@metadata$input)) {
    gi_obj@guideLFCs <- gi_obj@guideLFCs |>
      fill_gRNA_LFCs(gi_obj@metadata$input) |>
      collapse_replicates() |>
      compute_gene_main_effects()

    gi_obj@guideLFCs <- flatten_guide_lfcs(gi_obj@guideLFCs)
  }

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

      if (length(gi_obj@guideLFCs@data) > 0L) {
        gi_obj@guideLFCs@data[,, .r] <- make_symmetric(
          gi_obj@guideLFCs@data[,, .r]
        )
      }
    }

    if (identical(symmetric_analysis_method, "global_preaverage")) {
      gi_obj@guideGIs@data <- flatten_symmetric_pairs(
        .arr = gi_obj@guideGIs@data,
        pairs = gi_obj@screen_attr$unique_pairs
      )

      if (length(gi_obj@guideLFCs@data) > 0L) {
        gi_obj@guideLFCs@data <- flatten_symmetric_pairs(
          .arr = gi_obj@guideLFCs@data,
          pairs = gi_obj@screen_attr$unique_pairs
        )
      }

      gi_obj@guideGIs@space <- "gene_pair"
      gi_obj@guideLFCs@space <- "gene_pair"
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
