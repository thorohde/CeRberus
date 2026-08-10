make_screen_report_scores <- function() {
  genes <- paste0("G", seq_len(20L))
  input <- expand.grid(
    query_gene = genes,
    library_gene = genes,
    guide_pair = "g1",
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  input$GI <- seq_len(nrow(input))
  input
}

make_screen_for_report <- function(
  class = c("MultiplexScreen", "FixedPairScreen", "PosAgnMultiplexScreen"),
  checks = list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = FALSE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE
  ),
  failed_queries = character(),
  model_errors = list(),
  metadata = list(),
  limma_models = list()
) {
  base_screen <- GIScores(
    make_screen_report_scores(),
    block_layer = "guide_pair"
  )
  screen <- methods::as(base_screen, match.arg(class))
  selected_type <- if (methods::is(screen, "FixedPairScreen")) {
    "fixed_pair"
  } else {
    "multiplex"
  }

  screen@checks <- checks
  screen@errors$query_genes_not_usable <- failed_queries
  screen@errors$GI_computation_errors <- model_errors
  screen@metadata$requested_screen_type <- "auto"
  screen@metadata$inferred_screen_type <- selected_type
  screen@metadata$selected_screen_type <- selected_type
  screen@metadata <- utils::modifyList(screen@metadata, metadata)
  screen@limma_models <- limma_models

  if (methods::is(screen, "PosAgnMultiplexScreen")) {
    screen@metadata$symmetric_analysis_method <- "preaverage"
  }

  screen
}

add_screen_report_results <- function(screen) {
  if (methods::is(screen, "PosAgnMultiplexScreen")) {
    screen@symmGeneGIs <- data.table::data.table(
      gene_pair = c("A;B", "A;C", "B;C", "B;D", "C;D"),
      query_gene = c("A", "A", "B", "B", "C"),
      library_gene = c("B", "C", "C", "D", "D"),
      GI = c(1.2, -0.8, 0, NA_real_, Inf),
      GI_z = c(1, -1, 0, NA_real_, Inf),
      pval = c(0.01, 0.02, 0.5, NA_real_, 0.01),
      FDR = c(0.04, 0.05, 0.5, NA_real_, 0.01)
    )
  } else if (methods::is(screen, "FixedPairScreen")) {
    screen@geneGIs <- matrix(
      c(
        1.2, 0.01, 0.04,
        -0.8, 0.02, 0.05,
        0, 0.5, 0.5,
        NA_real_, NA_real_, NA_real_,
        Inf, 0.01, 0.01
      ),
      nrow = 5L,
      byrow = TRUE,
      dimnames = list(
        gene_pair = c("A;B", "A;C", "B;C", "B;D", "C;D"),
        variable = c("GI", "pval", "FDR")
      )
    )
  } else {
    screen@geneGIs <- array(
      c(
        1.2, -0.8, 0, NA_real_, Inf,
        0.01, 0.02, 0.5, NA_real_, 0.01,
        0.04, 0.05, 0.5, NA_real_, 0.01
      ),
      dim = c(1L, 5L, 3L),
      dimnames = list(
        query_gene = "A",
        library_gene = c("B", "C", "D", "E", "F"),
        variable = c("GI", "pval", "FDR")
      )
    )
  }

  screen@metadata$fdr_method <- "BH"
  screen
}

test_that("screen_report is part of the public package API", {
  expect_true("screen_report" %in% getNamespaceExports("CeRberus"))
})

test_that("screen_report returns a typed report", {
  screen <- make_screen_for_report(
    class = "MultiplexScreen",
    failed_queries = c("Q1", "Q2"),
    model_errors = list(Q1 = NULL, Q2 = simpleError("model failed"))
  )

  report <- screen_report(screen, print = FALSE)

  expect_type(report, "list")
  expect_named(
    report,
    c(
      "report_version",
      "screen",
      "dimensions",
      "model",
      "checks",
      "problems",
      "results"
    )
  )
  expect_identical(report$report_version, "1.0")
  expect_identical(report$screen$class, "MultiplexScreen")
  expect_identical(report$screen$interpreted_design, "multiplex")
  expect_identical(report$screen$requested_type, "auto")
  expect_identical(report$screen$inferred_type, "multiplex")
  expect_identical(report$screen$selected_type, "multiplex")
  expect_false(report$screen$position_agnostic)
  expect_null(report$screen$symmetric_analysis_method)

  expect_identical(report$dimensions$query_genes, 20L)
  expect_identical(report$dimensions$library_genes, 20L)
  expect_identical(report$dimensions$all_genes, 20L)
  expect_identical(report$dimensions$directional_pairs, 400L)
  expect_identical(report$dimensions$unordered_pairs, 210L)

  expect_identical(report$model$replicate_layers, "guide_pair")
  expect_identical(report$model$collapsed_layers, character())
  expect_identical(report$model$block_layer, "guide_pair")
  expect_identical(report$model$duplicate_correlation, numeric())
  expect_identical(report$model$fitted_models, 0L)

  expect_true(report$checks$gene_sets_equal)
  expect_true(report$checks$query_sufficient)
  expect_false(report$checks$library_sufficient)
  expect_true(report$checks$stable_library_size)
  expect_true(report$checks$sufficient_tests_per_query)

  expect_identical(report$problems$query_genes_missing_from_library, character())
  expect_identical(report$problems$library_genes_missing_from_query, character())
  expect_identical(report$problems$unusable_query_genes, c("Q1", "Q2"))
  expect_identical(unname(report$problems$model_errors), "model failed")

  expect_false(report$results$available)
  expect_null(report$results$fdr_method)
  expect_identical(report$results$fdr_threshold, 0.05)
  expect_null(report$results$tested_gene_pairs)
})

test_that("screen_report represents unavailable scalar values as NULL", {
  screen <- make_screen_for_report(
    class = "MultiplexScreen",
    checks = list()
  )
  screen@metadata$requested_screen_type <- NULL
  screen@metadata$inferred_screen_type <- NULL
  screen@metadata$selected_screen_type <- NULL
  methods::slot(screen@screen_attr, "n_query_genes", check = FALSE) <- NA_integer_
  methods::slot(screen@screen_attr, "n_lib_genes", check = FALSE) <- integer()
  methods::slot(screen@screen_attr, "n_all_genes", check = FALSE) <- NULL
  screen@guideGIs@block_layer <- character()

  report <- screen_report(screen, print = FALSE)

  expect_null(report$screen$requested_type)
  expect_null(report$screen$inferred_type)
  expect_null(report$screen$selected_type)
  expect_null(report$dimensions$query_genes)
  expect_null(report$dimensions$library_genes)
  expect_null(report$dimensions$all_genes)
  expect_null(report$model$block_layer)
  expect_true(all(purrr::map_lgl(report$checks, is.null)))
})

test_that("screen_report identifies every supported screen class", {
  expected_designs <- c(
    FixedPairScreen = "fixed-pair",
    MultiplexScreen = "multiplex",
    PosAgnMultiplexScreen = "position-agnostic multiplex"
  )

  purrr::iwalk(expected_designs, function(expected_design, class) {
    report <- screen_report(
      make_screen_for_report(class = class),
      print = FALSE
    )

    expect_identical(report$screen$class, class)
    expect_identical(report$screen$interpreted_design, expected_design)
    expect_identical(
      report$screen$position_agnostic,
      identical(class, "PosAgnMultiplexScreen")
    )

    if (identical(class, "PosAgnMultiplexScreen")) {
      expect_identical(report$screen$symmetric_analysis_method, "preaverage")
    } else {
      expect_null(report$screen$symmetric_analysis_method)
    }
  })
})

test_that("screen_report counts per-query and global models", {
  per_query_screen <- make_screen_for_report(
    class = "MultiplexScreen",
    limma_models = list(Q1 = list(), Q2 = NULL, Q3 = list())
  )
  global_fit <- limma::lmFit(matrix(c(1, 2, 3, 4), nrow = 2L))
  global_screen <- make_screen_for_report(
    class = "PosAgnMultiplexScreen",
    metadata = list(symmetric_analysis_method = "global_preaverage"),
    limma_models = global_fit
  )

  expect_identical(
    screen_report(per_query_screen, print = FALSE)$model$fitted_models,
    2L
  )
  expect_identical(
    screen_report(global_screen, print = FALSE)$model$fitted_models,
    1L
  )
})

test_that("screen_report summarizes results for every supported screen class", {
  for (class in c(
    "FixedPairScreen",
    "MultiplexScreen",
    "PosAgnMultiplexScreen"
  )) {
    report <- make_screen_for_report(class = class) |>
      add_screen_report_results() |>
      screen_report(print = FALSE)

    expect_true(report$results$available, info = class)
    expect_identical(report$results$fdr_method, "BH", info = class)
    expect_identical(report$results$fdr_threshold, 0.05, info = class)
    expected_tested_pairs <- if (identical(class, "MultiplexScreen")) {
      4L
    } else {
      5L
    }
    expect_identical(
      report$results$tested_gene_pairs,
      expected_tested_pairs,
      info = class
    )
    expect_identical(report$results$finite_results, 3L, info = class)
    expect_identical(report$results$significant_results, 3L, info = class)
    expect_identical(
      report$results$positive_significant_results,
      1L,
      info = class
    )
    expect_identical(
      report$results$negative_significant_results,
      1L,
      info = class
    )
  }
})

test_that("screen_report prints a readable report", {
  screen <- make_screen_for_report(
    class = "MultiplexScreen",
    failed_queries = paste0("Q", seq_len(10L))
  )

  output <- capture.output(screen_report(screen, print = TRUE, width = 12))
  text <- paste(output, collapse = "\n")

  expect_equal(output[[2L]], paste(rep("=", 12L), collapse = ""))
  expect_match(text, "CeRberus screen report", fixed = TRUE)
  expect_match(text, "OVERVIEW", fixed = TRUE)
  expect_match(text, "MODEL", fixed = TRUE)
  expect_match(text, "CHECKS", fixed = TRUE)
  expect_match(text, "PROBLEMS", fixed = TRUE)
  expect_match(text, "RESULTS", fixed = TRUE)
  expect_match(text, "query_sufficient: OK", fixed = TRUE)
  expect_match(text, "library_sufficient: PROBLEM", fixed = TRUE)
  expect_match(text, "+2 more", fixed = TRUE)
})

test_that("screen_report writes YAML and suppresses printing by default", {
  screen <- make_screen_for_report(class = "MultiplexScreen")

  for (extension in c(".yaml", ".yml")) {
    path <- file.path(tempdir(), paste0("nested-report", extension))
    unlink(path)

    output <- capture.output(
      report <- screen_report(screen, file = path)
    )
    written <- yaml::read_yaml(path)

    expect_length(output, 0L)
    expect_true(file.exists(path))
    expect_identical(report$screen$class, "MultiplexScreen")
    expect_identical(written$report_version, "1.0")
    expect_identical(written$screen$class, "MultiplexScreen")
    expect_false(written$screen$position_agnostic)
    expect_identical(written$dimensions$query_genes, 20L)
    expect_false(written$checks$library_sufficient)
  }
})

test_that("screen_report creates the report parent directory", {
  screen <- make_screen_for_report(class = "MultiplexScreen")
  directory <- tempfile("screen-report-directory-")
  path <- file.path(directory, "nested", "screen_report.yaml")

  screen_report(screen, file = path, print = FALSE)

  expect_true(file.exists(path))
})

test_that("screen_report validates output controls", {
  screen <- make_screen_for_report(class = "MultiplexScreen")

  invalid_files <- list(NA_character_, character(), c("a.yaml", "b.yaml"), 1)
  purrr::walk(invalid_files, function(value) {
    expect_error(
      screen_report(screen, file = value, print = FALSE),
      "file must be NULL or a single path string"
    )
  })
  expect_error(
    screen_report(screen, file = "report.txt", print = FALSE),
    "file must use a .yaml or .yml extension"
  )

  invalid_logical_values <- list(NA, logical(), c(TRUE, FALSE), 1, "TRUE")
  purrr::walk(invalid_logical_values, function(value) {
    expect_error(
      screen_report(screen, print = value),
      "print must be TRUE or FALSE"
    )
  })

  invalid_widths <- list(NA_real_, Inf, 0, -1, numeric(), c(40, 80), "80")
  purrr::walk(invalid_widths, function(value) {
    expect_error(
      screen_report(screen, print = FALSE, width = value),
      "width must be one finite positive numeric value"
    )
  })
})
