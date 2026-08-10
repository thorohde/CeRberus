#####

.screen_report_scalar <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(NULL)
  }

  x[[1L]]
}

.screen_report_model_count <- function(gi_obj) {
  if (inherits(gi_obj@limma_models, "MArrayLM")) {
    return(1L)
  }

  as.integer(sum(!purrr::map_lgl(gi_obj@limma_models, is.null)))
}

.screen_report_error_messages <- function(errors) {
  errors <- purrr::compact(errors)

  if (length(errors) == 0L) {
    return(character())
  }

  purrr::map_chr(errors, function(error) {
    if (inherits(error, "condition")) {
      conditionMessage(error)
    } else {
      as.character(error)[[1L]]
    }
  })
}

.screen_report_results <- function(gi_obj, fdr_threshold = 0.05) {
  results_available <- length(gi_obj@geneGIs) > 0L ||
    (
      methods::is(gi_obj, "PosAgnMultiplexScreen") &&
        nrow(gi_obj@symmGeneGIs) > 0L
    )

  output <- list(
    available = results_available,
    fdr_method = .screen_report_scalar(gi_obj@metadata$fdr_method),
    fdr_threshold = fdr_threshold,
    tested_gene_pairs = NULL,
    finite_results = NULL,
    significant_results = NULL,
    positive_significant_results = NULL,
    negative_significant_results = NULL
  )

  if (!results_available) {
    return(output)
  }

  results <- gi_df(gi_obj)
  required_columns <- c("GI", "FDR")

  if (!all(required_columns %in% colnames(results))) {
    return(output)
  }

  finite_gi <- is.finite(results$GI)
  significant <- is.finite(results$FDR) & results$FDR <= fdr_threshold

  output$tested_gene_pairs <- as.integer(nrow(results))
  output$finite_results <- as.integer(sum(finite_gi))
  output$significant_results <- as.integer(sum(significant))
  output$positive_significant_results <- as.integer(sum(
    significant & finite_gi & results$GI > 0
  ))
  output$negative_significant_results <- as.integer(sum(
    significant & finite_gi & results$GI < 0
  ))

  output
}

.build_screen_report <- function(gi_obj) {
  screen_class <- class(gi_obj)[[1L]]
  screen_types <- c(
    FixedPairScreen = "fixed-pair",
    MultiplexScreen = "multiplex",
    PosAgnMultiplexScreen = "position-agnostic multiplex"
  )
  interpreted_design <- unname(screen_types[screen_class])
  if (is.na(interpreted_design)) {
    interpreted_design <- screen_class
  }

  screen_attributes <- gi_obj@screen_attr
  metadata <- gi_obj@metadata
  errors <- gi_obj@errors
  failed_queries <- errors$query_genes_not_usable
  if (is.null(failed_queries)) {
    failed_queries <- character()
  }

  symmetric_analysis_method <- NULL
  if (methods::is(gi_obj, "PosAgnMultiplexScreen")) {
    symmetric_analysis_method <- get_symmetric_analysis_method(gi_obj)
  }

  list(
    report_version = "1.0",
    screen = list(
      class = screen_class,
      interpreted_design = interpreted_design,
      requested_type = .screen_report_scalar(metadata$requested_screen_type),
      inferred_type = .screen_report_scalar(metadata$inferred_screen_type),
      selected_type = .screen_report_scalar(metadata$selected_screen_type),
      position_agnostic = methods::is(gi_obj, "PosAgnMultiplexScreen"),
      symmetric_analysis_method = symmetric_analysis_method
    ),
    dimensions = list(
      query_genes = .screen_report_scalar(screen_attributes@n_query_genes),
      library_genes = .screen_report_scalar(screen_attributes@n_lib_genes),
      all_genes = .screen_report_scalar(screen_attributes@n_all_genes),
      directional_pairs = as.integer(length(screen_attributes@all_pairs)),
      unordered_pairs = as.integer(length(screen_attributes@unique_pairs))
    ),
    model = list(
      replicate_layers = as.character(gi_obj@guideGIs@replicates),
      collapsed_layers = as.character(gi_obj@guideGIs@collapse),
      block_layer = .screen_report_scalar(gi_obj@guideGIs@block_layer),
      duplicate_correlation = as.numeric(gi_obj@dupCorrelation),
      fitted_models = .screen_report_model_count(gi_obj)
    ),
    checks = list(
      gene_sets_equal = .screen_report_scalar(gi_obj@checks$gene_sets_equal),
      query_sufficient = .screen_report_scalar(gi_obj@checks$query_sufficient),
      library_sufficient = .screen_report_scalar(gi_obj@checks$library_sufficient),
      stable_library_size = .screen_report_scalar(
        gi_obj@checks$stable_library_size
      ),
      sufficient_tests_per_query = .screen_report_scalar(
        gi_obj@checks$sufficient_tests_per_query
      )
    ),
    problems = list(
      query_genes_missing_from_library = as.character(
        screen_attributes@query_genes_not_in_lib
      ),
      library_genes_missing_from_query = as.character(
        screen_attributes@library_genes_not_in_query
      ),
      unusable_query_genes = as.character(failed_queries),
      model_errors = .screen_report_error_messages(errors$GI_computation_errors)
    ),
    results = .screen_report_results(gi_obj)
  )
}

.print_screen_report <- function(report, width = 80) {
  line <- function(character = "-") {
    paste(rep(character, width), collapse = "")
  }
  display_value <- function(value) {
    if (is.null(value) || length(value) == 0L) {
      return("not available")
    }
    paste(value, collapse = ", ")
  }
  display_check <- function(value) {
    if (is.null(value)) {
      return("not checked")
    }
    if (isTRUE(value)) {
      return("OK")
    }
    "PROBLEM"
  }
  display_items <- function(items, n = 8L) {
    if (length(items) == 0L) {
      return("none")
    }

    suffix <- if (length(items) > n) {
      paste0(" ... +", length(items) - n, " more")
    } else {
      ""
    }
    paste0(paste(utils::head(items, n), collapse = ", "), suffix)
  }

  sections <- list(
    overview = c(
      paste0("Screen class: ", display_value(report$screen$class)),
      paste0(
        "Interpreted design: ",
        display_value(report$screen$interpreted_design)
      ),
      paste0("Query genes: ", display_value(report$dimensions$query_genes)),
      paste0("Library genes: ", display_value(report$dimensions$library_genes)),
      paste0("All genes: ", display_value(report$dimensions$all_genes)),
      paste0(
        "Directional pairs: ",
        display_value(report$dimensions$directional_pairs)
      ),
      paste0(
        "Unordered pairs: ",
        display_value(report$dimensions$unordered_pairs)
      )
    ),
    model = c(
      paste0(
        "Replicate layers: ",
        display_value(report$model$replicate_layers)
      ),
      paste0(
        "Collapsed layers: ",
        display_value(report$model$collapsed_layers)
      ),
      paste0("Block layer: ", display_value(report$model$block_layer)),
      paste0(
        "Duplicate correlation: ",
        display_value(report$model$duplicate_correlation)
      ),
      paste0("Fitted models: ", display_value(report$model$fitted_models))
    ),
    checks = purrr::imap_chr(report$checks, function(value, name) {
      paste0(name, ": ", display_check(value))
    }),
    problems = c(
      paste0(
        "Query genes missing from library: ",
        display_items(report$problems$query_genes_missing_from_library)
      ),
      paste0(
        "Library genes missing from query: ",
        display_items(report$problems$library_genes_missing_from_query)
      ),
      paste0(
        "Unusable query genes: ",
        display_items(report$problems$unusable_query_genes)
      ),
      paste0(
        "Model errors: ",
        display_items(report$problems$model_errors)
      )
    ),
    results = purrr::imap_chr(report$results, function(value, name) {
      paste0(name, ": ", display_value(value))
    })
  )

  cat("\n", line("="), "\n", sep = "")
  cat("CeRberus screen report\n")
  cat(line("="), "\n", sep = "")

  purrr::iwalk(sections, function(section, name) {
    cat("\n", toupper(name), "\n", sep = "")
    cat(line("-"), "\n", sep = "")
    cat(paste0("- ", section, collapse = "\n"), "\n", sep = "")
  })

  invisible(report)
}

.write_screen_report_yaml <- function(report, file) {
  stopifnot(
    "file must be a single path string." = is.character(file) &&
      length(file) == 1L &&
      !is.na(file) &&
      nzchar(file),
    "file must use a .yaml or .yml extension." = tolower(tools::file_ext(file)) %in%
      c("yaml", "yml")
  )

  parent_directory <- dirname(file)
  if (!dir.exists(parent_directory)) {
    dir.create(parent_directory, recursive = TRUE)
  }

  yaml::write_yaml(report, file = file)
  invisible(file)
}

#' @describeIn screen_report Summarize a `ScreenBase` object.
#'
#' @param file Optional path to a `.yaml` or `.yml` report file.
#' @param print Logical scalar. Print the report to the console when `TRUE`.
#' @param width Positive numeric scalar controlling the printed separator width.
setMethod(
  "screen_report",
  signature = signature(gi_obj = "ScreenBase"),
  function(gi_obj, file = NULL, print = is.null(file), width = 80) {
    stopifnot(
      "file must be NULL or a single path string." = is.null(file) ||
        (
          is.character(file) &&
            length(file) == 1L &&
            !is.na(file) &&
            nzchar(file)
        ),
      "print must be TRUE or FALSE." = is.logical(print) &&
        length(print) == 1L &&
        !is.na(print),
      "width must be one finite positive numeric value." = is.numeric(width) &&
        length(width) == 1L &&
        is.finite(width) &&
        width > 0
    )

    report <- .build_screen_report(gi_obj)

    if (!is.null(file)) {
      .write_screen_report_yaml(report, file)
    }

    if (isTRUE(print)) {
      .print_screen_report(report, width = width)
    }

    report
  }
)

#####