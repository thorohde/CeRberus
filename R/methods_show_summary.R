#' Display and summarize CeRberus objects
#'
#' @description
#' `show()` provides a compact display of CeRberus S4 objects without printing
#' their potentially large arrays, fitted models, or metadata. `summary()`
#' returns a structured overview with data dimensions, completeness statistics,
#' and class-specific configuration.
#'
#' Screen subclasses inherit both methods from [`ScreenBase-class`]. Their
#' concrete class and result representation are reflected in the output.
#'
#' @param object A CeRberus S4 object.
#' @param x A summary returned by one of the `summary()` methods.
#' @param ... Additional arguments, currently unused.
#'
#' @return `show()` returns `object` invisibly. `summary()` returns a classed
#'   list. The summary print method returns `x` invisibly.
#'
#' @name show-summary
NULL

.display_dimensions <- function(x) {
  dimensions <- dim(x)
  if (is.null(dimensions) || length(dimensions) == 0L) {
    return("empty")
  }

  paste(dimensions, collapse = " x ")
}

.display_items <- function(x, empty = "none") {
  if (is.null(x) || length(x) == 0L) {
    return(empty)
  }

  paste(x, collapse = ", ")
}

.numeric_summary <- function(x) {
  finite_values <- as.numeric(x)[is.finite(x)]

  list(
    values = as.integer(length(x)),
    missing = as.integer(sum(is.na(x))),
    finite = as.integer(sum(is.finite(x))),
    non_finite = as.integer(sum(!is.na(x) & !is.finite(x))),
    minimum = if (length(finite_values) > 0L) min(finite_values) else NULL,
    median = if (length(finite_values) > 0L) {
      stats::median(finite_values)
    } else {
      NULL
    },
    mean = if (length(finite_values) > 0L) mean(finite_values) else NULL,
    maximum = if (length(finite_values) > 0L) max(finite_values) else NULL
  )
}

.array_summary <- function(x) {
  c(
    list(
      dimensions = as.integer(dim(x)),
      dimension_names = names(dimnames(x))
    ),
    .numeric_summary(x)
  )
}

.effect_summary <- function(x) {
  list(
    total = as.integer(length(x)),
    missing = as.integer(sum(is.na(x))),
    finite = as.integer(sum(is.finite(x)))
  )
}

.new_cerberus_summary <- function(x, class) {
  structure(x, class = c(class, "CeRberusSummary", "list"))
}

.screen_model_count <- function(models) {
  if (inherits(models, "MArrayLM")) {
    return(1L)
  }
  if (!is.list(models)) {
    return(as.integer(length(models) > 0L))
  }

  as.integer(sum(!purrr::map_lgl(models, is.null)))
}

.screen_display_slots <- function(object) {
  if (methods::is(object, "PosAgnMultiplexScreen")) {
    return(list(
      guide_gis = object@aggregatedGuideGIs,
      models = object@aggregatedLimmaModels,
      failed_queries = object@errors$aggregated_query_genes_not_usable,
      model_errors = object@errors$aggregated_GI_computation_errors
    ))
  }

  list(
    guide_gis = object@guideGIs,
    models = object@limma_models,
    failed_queries = object@errors$query_genes_not_usable,
    model_errors = object@errors$GI_computation_errors
  )
}

.condition_count <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(0L)
  }
  if (inherits(x, "condition") || is.character(x)) {
    return(as.integer(length(x)))
  }
  if (is.list(x)) {
    return(as.integer(sum(purrr::map_int(x, .condition_count))))
  }

  0L
}

.screen_result_description <- function(object) {
  if (methods::is(object, "PosAgnMultiplexScreen")) {
    if (nrow(object@symmGeneGIs) == 0L) {
      return("not available")
    }
    return(paste(nrow(object@symmGeneGIs), "symmetric gene pairs"))
  }

  if (length(object@geneGIs) == 0L) {
    return("not available")
  }

  paste0("available (", .display_dimensions(object@geneGIs), ")")
}

.summary_label <- function(x) {
  tools::toTitleCase(gsub("_", " ", x, fixed = TRUE))
}

.summary_value <- function(x, max_items = 8L) {
  if (is.null(x)) {
    return("not available")
  }
  if (length(x) == 0L) {
    return("none")
  }

  shown <- utils::head(x, max_items)
  if (is.numeric(shown)) {
    shown <- signif(shown, 5L)
  }
  suffix <- if (length(x) > max_items) {
    paste0(" ... +", length(x) - max_items, " more")
  } else {
    ""
  }

  paste0(paste(shown, collapse = ", "), suffix)
}

.print_summary_entries <- function(x, indent = 0L) {
  indentation <- paste(rep(" ", indent), collapse = "")

  purrr::iwalk(x, function(value, name) {
    label <- .summary_label(name)
    if (is.list(value) && !inherits(value, "data.frame")) {
      cat(indentation, label, "\n", sep = "")
      .print_summary_entries(value, indent = indent + 2L)
    } else {
      cat(
        indentation,
        label,
        ": ",
        .summary_value(value),
        "\n",
        sep = ""
      )
    }
  })
}

#' @rdname show-summary
#' @aliases show,gRNA_LFC-method
setMethod("show", "gRNA_LFC", function(object) {
  cat("<gRNA_LFC>\n")
  cat("  Data: ", .display_dimensions(object@data), "\n", sep = "")
  cat("  Space: ", .display_items(object@space), "\n", sep = "")
  cat("  Replicates: ", .display_items(object@replicates), "\n", sep = "")
  cat("  Collapsed: ", .display_items(object@collapse), "\n", sep = "")
  cat(
    "  Main effects: ",
    length(object@query_main_effects),
    " query, ",
    length(object@library_main_effects),
    " library\n",
    sep = ""
  )
  invisible(object)
})

#' @rdname show-summary
#' @aliases summary,gRNA_LFC-method
setMethod("summary", "gRNA_LFC", function(object, ...) {
  .new_cerberus_summary(
    list(
      class = class(object)[[1L]],
      data = .array_summary(object@data),
      space = object@space,
      replicates = object@replicates,
      collapsed_layers = object@collapse,
      main_effects = list(
        query = .effect_summary(object@query_main_effects),
        library = .effect_summary(object@library_main_effects)
      )
    ),
    "summary_gRNA_LFC"
  )
})

#' @rdname show-summary
#' @aliases show,gRNA_GI-method
setMethod("show", "gRNA_GI", function(object) {
  cat("<gRNA_GI>\n")
  cat("  Data: ", .display_dimensions(object@data), "\n", sep = "")
  cat("  Space: ", .display_items(object@space), "\n", sep = "")
  cat("  Replicates: ", .display_items(object@replicates), "\n", sep = "")
  cat("  Collapsed: ", .display_items(object@collapse), "\n", sep = "")
  cat(
    "  Blocking: ",
    if (isTRUE(object@use_blocks)) "active" else "inactive",
    if (length(object@block_layer) > 0L) {
      paste0(" (", object@block_layer[[1L]], ")")
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  invisible(object)
})

#' @rdname show-summary
#' @aliases summary,gRNA_GI-method
setMethod("summary", "gRNA_GI", function(object, ...) {
  .new_cerberus_summary(
    list(
      class = class(object)[[1L]],
      data = .array_summary(object@data),
      space = object@space,
      replicates = object@replicates,
      collapsed_layers = object@collapse,
      blocking = list(
        active = isTRUE(object@use_blocks),
        layer = if (length(object@block_layer) > 0L) {
          object@block_layer[[1L]]
        } else {
          NULL
        },
        assignments = as.integer(length(object@blocks)),
        descriptions = as.integer(length(object@block_description))
      )
    ),
    "summary_gRNA_GI"
  )
})

#' @rdname show-summary
#' @aliases show,ScreenDesign-method
setMethod("show", "ScreenDesign", function(object) {
  cat("<ScreenDesign>\n")
  cat(
    "  Genes: ",
    object@n_query_genes,
    " query, ",
    object@n_lib_genes,
    " library, ",
    object@n_all_genes,
    " total\n",
    sep = ""
  )
  cat(
    "  Pairs: ",
    length(object@all_pairs),
    " directional, ",
    length(object@unique_pairs),
    " unordered\n",
    sep = ""
  )
  cat(
    "  Axis-specific genes: ",
    length(object@query_genes_not_in_lib),
    " query, ",
    length(object@library_genes_not_in_query),
    " library\n",
    sep = ""
  )
  invisible(object)
})

#' @rdname show-summary
#' @aliases summary,ScreenDesign-method
setMethod("summary", "ScreenDesign", function(object, ...) {
  .new_cerberus_summary(
    list(
      class = class(object)[[1L]],
      genes = list(
        query = object@n_query_genes,
        library = object@n_lib_genes,
        total = object@n_all_genes,
        query_only = as.integer(length(object@query_genes_not_in_lib)),
        library_only = as.integer(length(object@library_genes_not_in_query))
      ),
      pairs = list(
        directional = as.integer(length(object@all_pairs)),
        unordered = as.integer(length(object@unique_pairs))
      ),
      observations_per_query = .numeric_summary(
        object@observations_per_query
      ),
      contrasts_present = !is.null(object@contrasts) &&
        length(object@contrasts) > 0L
    ),
    "summary_ScreenDesign"
  )
})

#' @rdname show-summary
#' @aliases show,ScreenBase-method
setMethod("show", "ScreenBase", function(object) {
  screen_types <- c(
    FixedPairScreen = "fixed-pair",
    MultiplexScreen = "multiplex",
    PosAgnMultiplexScreen = "position-agnostic multiplex",
    ScreenBase = "base"
  )
  screen_class <- class(object)[[1L]]
  design <- unname(screen_types[screen_class])
  if (is.na(design)) {
    design <- screen_class
  }

  slots <- .screen_display_slots(object)
  failed_queries <- slots$failed_queries
  if (is.null(failed_queries)) {
    failed_queries <- character()
  }

  cat("<", screen_class, ">\n", sep = "")
  cat("  Design: ", design, "\n", sep = "")
  cat(
    "  Genes: ",
    object@screen_attr@n_query_genes,
    " query, ",
    object@screen_attr@n_lib_genes,
    " library, ",
    object@screen_attr@n_all_genes,
    " total\n",
    sep = ""
  )
  cat(
    "  Pairs: ",
    length(object@screen_attr@all_pairs),
    " directional, ",
    length(object@screen_attr@unique_pairs),
    " unordered\n",
    sep = ""
  )
  cat(
    "  Guide GI data: ",
    .display_dimensions(slots$guide_gis@data),
    "\n",
    sep = ""
  )
  cat("  Models: ", .screen_model_count(slots$models), " fitted\n", sep = "")
  cat("  Results: ", .screen_result_description(object), "\n", sep = "")
  cat(
    "  Problems: ",
    length(failed_queries),
    " unusable queries, ",
    .condition_count(slots$model_errors),
    " model errors\n",
    sep = ""
  )
  invisible(object)
})

#' @rdname show-summary
#' @aliases summary,ScreenBase-method
setMethod("summary", "ScreenBase", function(object, ...) {
  .new_cerberus_summary(
    .build_screen_report(object),
    "summary_CeRberusScreen"
  )
})

#' @rdname show-summary
#' @export
print.CeRberusSummary <- function(x, ...) {
  cat("<", class(x)[[1L]], ">\n", sep = "")
  .print_summary_entries(unclass(x))
  invisible(x)
}
