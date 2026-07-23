#' @describeIn collapse_replicates Collapse replicate dimensions in a
#'   [`gRNA_GI-class`] object by averaging over layers listed in the `collapse`
#'   slot.
#' @aliases collapse_replicates,gRNA_GI-method

#####

setMethod(
  "collapse_replicates",
  signature = signature(.x = "gRNA_GI"),
  function(.x) {
    if (length(.x@collapse) == 0) {
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

    collapsed <- apply(
      X = .x@data,
      MARGIN = keep_margin,
      FUN = mean,
      na.rm = TRUE
    )

    # Rebuild array explicitly so dimension order and dimnames match keep_dims.
    keep_dimnames <- dimnames(.x@data)[keep_margin]
    names(keep_dimnames) <- keep_dims

    collapsed <- array(
      data = collapsed,
      dim = map_int(keep_dimnames, length),
      dimnames = keep_dimnames
    )

    .x@data <- collapsed

    .x@replicates <- setdiff(.x@replicates, .x@collapse)

    return(.x)
  }
)

#####

setMethod(
  "flatten_guideGIs",
  signature = signature(.x = "gRNA_GI"),
  function(.x) {
    .x@use_blocks <- .x@block_layer != "" &
      length(.x@block_layer) > 0 &
      length(setdiff(.x@replicates, .x@block_layer)) != 0

    .f <- as.formula(paste0(
      paste0(.x@space, collapse = " ~ "),
      " ~ ",
      paste0(.x@replicates, collapse = " + ")
    ))

    .x@data <- .x@data |>
      flatten_array(c(.x@space, .x@replicates)) |>
      acast(.f)

    .x@block_description <- dimnames(.x@data)[[length(.x@space) + 1]]

    if (length(.x@blocks) != 1 && !identical(.x@blocks, "none")) {
      .x@blocks <- c(
        guide_pair = "(g\\d+)",
        bio_rep = "(b\\d+)",
        tech_rep = "(t\\d+)"
      )[[.x@block_layer]]

      .x@blocks <- as.character(factor(str_match(
        .x@block_description,
        .x@blocks
      )[, 2]))
    }
    return(.x)
  }
)

#####

setMethod(
  "compute_dupCorrelation",
  signature = signature(.x = "gRNA_GI"),
  function(.x) {
    if (length(.x@space) == 1) {
      if (.x@use_blocks) {
        output <- suppressWarnings(
          limma::duplicateCorrelation(
            object = .x@data,
            block = .x@blocks
          )
        )
      } else {
        output <- suppressWarnings(
          limma::duplicateCorrelation(object = .x@data, ndups = 1)
        )
      }
      output <- output$consensus.correlation
    }

    if (length(.x@space) == 2) {
      output <- set_names(rownames(.x@data)) |>
        map(\(.g) .x@data[.g, , ])

      if (.x@use_blocks) {
        output <- output |>
          map(\(.g) {
            suppressWarnings(
              limma::duplicateCorrelation(
                object = .g,
                block = .x@blocks
              )
            )
          })
      } else {
        output <- output |>
          map(\(.g) {
            suppressWarnings(
              limma::duplicateCorrelation(
                object = .g,
                ndups = 1 # required
              )
            )
          })
      }
      output <- output |> map_dbl("consensus.correlation")
    }
    return(output)
  }
)

#####
