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

make_combined_report_screens <- function() {
  duplicate_correlation_data <- data.table::data.table(
    config = c("configuration_a", "configuration_b"),
    dcor = c(0.1, 0.2),
    kept = c("selected", "")
  )

  screens <- list(
    configuration_a = make_screen_for_report(
      class = "MultiplexScreen",
      limma_models = list(model = list())
    ),
    configuration_b = make_screen_for_report(class = "FixedPairScreen")
  )

  purrr::map(screens, function(screen) {
    screen@metadata$dupcor_data <- duplicate_correlation_data
    screen
  })
}

test_that("screen_report exposes a stable schema and scalar representations", {
  screen <- make_screen_for_report(
    class = "MultiplexScreen",
    failed_queries = c("Q1", "Q2"),
    model_errors = list(Q1 = NULL, Q2 = simpleError("model failed"))
  )
  screen@dupCorrelation <- c(Q1 = 0.1, Q2 = 0.2, Q3 = NA_real_)

  report <- screen_report(screen, print = FALSE)

  expect_equal(
    list(
      exported = "screen_report" %in% getNamespaceExports("CeRberus"),
      report = report
    ),
    list(
      exported = TRUE,
      report = list(
        report_version = "1.1",
        screen = list(
          class = "MultiplexScreen",
          interpreted_design = "multiplex",
          requested_type = "auto",
          inferred_type = "multiplex",
          selected_type = "multiplex",
          position_agnostic = FALSE,
          symmetric_analysis_method = NULL
        ),
        dimensions = list(
          query_genes = 20L,
          library_genes = 20L,
          all_genes = 20L,
          directional_pairs = 400L,
          unordered_pairs = 210L
        ),
        model = list(
          replicate_layers = "guide_pair",
          collapsed_layers = character(),
          block_layer = "guide_pair",
          duplicate_correlation = 0.15,
          fitted_models = 0L
        ),
        checks = list(
          gene_sets_equal = TRUE,
          query_sufficient = TRUE,
          library_sufficient = FALSE,
          stable_library_size = TRUE,
          sufficient_tests_per_query = TRUE
        ),
        problems = list(
          query_genes_missing_from_library = character(),
          library_genes_missing_from_query = character(),
          unusable_query_genes = c("Q1", "Q2"),
          model_errors = c(Q2 = "model failed")
        ),
        results = list(
          available = FALSE,
          genomic_inflation = NULL,
          fdr_method = NULL,
          fdr_threshold = 0.05,
          tested_gene_pairs = NULL,
          finite_results = NULL,
          significant_results = NULL,
          positive_significant_results = NULL,
          negative_significant_results = NULL
        )
      )
    )
  )

  unavailable_screen <- make_screen_for_report(
    class = "MultiplexScreen",
    checks = list()
  )
  unavailable_screen@metadata$requested_screen_type <- NULL
  unavailable_screen@metadata$inferred_screen_type <- NULL
  unavailable_screen@metadata$selected_screen_type <- NULL
  methods::slot(
    unavailable_screen@screen_attr,
    "n_query_genes",
    check = FALSE
  ) <- NA_integer_
  methods::slot(
    unavailable_screen@screen_attr,
    "n_lib_genes",
    check = FALSE
  ) <- integer()
  methods::slot(
    unavailable_screen@screen_attr,
    "n_all_genes",
    check = FALSE
  ) <- NULL
  unavailable_screen@guideGIs@block_layer <- character()

  unavailable_report <- screen_report(unavailable_screen, print = FALSE)

  expect_equal(
    list(
      screen = unavailable_report$screen[c(
        "requested_type",
        "inferred_type",
        "selected_type"
      )],
      dimensions = unavailable_report$dimensions[c(
        "query_genes",
        "library_genes",
        "all_genes"
      )],
      block_layer = unavailable_report$model["block_layer"],
      checks = unavailable_report$checks
    ),
    list(
      screen = list(
        requested_type = NULL,
        inferred_type = NULL,
        selected_type = NULL
      ),
      dimensions = list(
        query_genes = NULL,
        library_genes = NULL,
        all_genes = NULL
      ),
      block_layer = list(block_layer = NULL),
      checks = list(
        gene_sets_equal = NULL,
        query_sufficient = NULL,
        library_sufficient = NULL,
        stable_library_size = NULL,
        sufficient_tests_per_query = NULL
      )
    )
  )
})

test_that("screen_report identifies screen classes and counts model strategies", {
  global_fit <- limma::lmFit(matrix(c(1, 2, 3, 4), nrow = 2L))
  global_screen <- make_screen_for_report(
    class = "PosAgnMultiplexScreen",
    limma_models = global_fit
  )
  global_screen@metadata$symmetric_analysis_method <- "global_preaverage"
  cases <- list(
    FixedPairScreen = list(
      screen = make_screen_for_report(class = "FixedPairScreen"),
      design = "fixed-pair",
      position_agnostic = FALSE,
      symmetric_analysis_method = NULL,
      fitted_models = 0L
    ),
    MultiplexScreen = list(
      screen = make_screen_for_report(
        class = "MultiplexScreen",
        limma_models = list(Q1 = list(), Q2 = NULL, Q3 = list())
      ),
      design = "multiplex",
      position_agnostic = FALSE,
      symmetric_analysis_method = NULL,
      fitted_models = 2L
    ),
    PosAgnMultiplexScreen = list(
      screen = global_screen,
      design = "position-agnostic multiplex",
      position_agnostic = TRUE,
      symmetric_analysis_method = "global_preaverage",
      fitted_models = 1L
    )
  )

  purrr::iwalk(cases, function(case, class) {
    report <- screen_report(case$screen, print = FALSE)

    expect_equal(
      list(
        class = report$screen$class,
        design = report$screen$interpreted_design,
        position_agnostic = report$screen$position_agnostic,
        symmetric_analysis_method = report$screen$symmetric_analysis_method,
        fitted_models = report$model$fitted_models
      ),
      list(
        class = class,
        design = case$design,
        position_agnostic = case$position_agnostic,
        symmetric_analysis_method = case$symmetric_analysis_method,
        fitted_models = case$fitted_models
      ),
      info = class
    )
  })
})

test_that("screen_report summarizes results for every supported screen class", {
  expected_genomic_inflation <- stats::median(stats::qchisq(
    1 - c(0.01, 0.02, 0.5, 0.01),
    df = 1
  )) / stats::qchisq(0.5, df = 1)

  for (class in c(
    "FixedPairScreen",
    "MultiplexScreen",
    "PosAgnMultiplexScreen"
  )) {
    report <- make_screen_for_report(class = class) |>
      add_screen_report_results() |>
      screen_report(print = FALSE)

    expected_tested_pairs <- if (identical(class, "MultiplexScreen")) {
      4L
    } else {
      5L
    }
    expect_equal(
      report$results,
      list(
        available = TRUE,
        genomic_inflation = expected_genomic_inflation,
        fdr_method = "BH",
        fdr_threshold = 0.05,
        tested_gene_pairs = expected_tested_pairs,
        finite_results = 3L,
        significant_results = 3L,
        positive_significant_results = 1L,
        negative_significant_results = 1L
      ),
      info = class
    )
  }
})

test_that("combined screen report preserves structure through YAML", {
  screens <- make_combined_report_screens()
  screens$configuration_a <- add_screen_report_results(
    screens$configuration_a
  )

  report <- build_combined_screen_report(screens)

  path <- tempfile(fileext = ".yaml")
  write_screen_report_yaml(report, path)
  written <- yaml::read_yaml(path)

  expected_genomic_inflation <- stats::median(stats::qchisq(
    1 - c(0.01, 0.02, 0.5, 0.01),
    df = 1
  )) / stats::qchisq(0.5, df = 1)
  report_contract <- function(value) {
    list(
      names = names(value),
      report_version = value$report_version,
      selection = value$selection,
      configuration_names = names(value$configurations),
      computed = list(
        configuration = value$configurations$configuration_a$screen$configuration,
        class = value$configurations$configuration_a$screen$class,
        genomic_inflation = value$configurations$configuration_a$results$genomic_inflation,
        report_version = value$configurations$configuration_a$report_version
      ),
      not_computed = value$configurations$configuration_b
    )
  }
  expected <- list(
    names = c("report_version", "selection", "configurations"),
    report_version = "1.1",
    selection = list(
      selected_configuration = "configuration_a",
      evaluated_configurations = list(
        configuration_a = list(duplicate_correlation = 0.1, selected = TRUE),
        configuration_b = list(duplicate_correlation = 0.2, selected = FALSE)
      )
    ),
    configuration_names = c("configuration_a", "configuration_b"),
    computed = list(
      configuration = "configuration_a",
      class = "MultiplexScreen",
      genomic_inflation = expected_genomic_inflation,
      report_version = NULL
    ),
    not_computed = list(
      status = "not_computed",
      screen = list(configuration = "configuration_b"),
      model = list(duplicate_correlation = NULL)
    )
  )

  expect_equal(report_contract(report), expected)
  expect_equal(report_contract(written), expected)
})

test_that("combined screen report validates its screen-object list", {
  screen <- make_screen_for_report(class = "MultiplexScreen")
  cases <- list(
    empty = list(
      input = list(),
      error = "screen_objects must be a non-empty list"
    ),
    unnamed = list(
      input = list(screen),
      error = "screen_objects must be named"
    ),
    empty_name = list(
      input = stats::setNames(list(screen), ""),
      error = "screen_objects must be named"
    ),
    duplicate_names = list(
      input = structure(
        list(screen, screen),
        names = c("configuration", "configuration")
      ),
      error = "screen object names must be unique"
    ),
    invalid_entry = list(
      input = list(configuration = list()),
      error = "Every entry must inherit from ScreenBase"
    )
  )

  purrr::iwalk(cases, function(case, name) {
    expect_error(
      build_combined_screen_report(case$input),
      case$error,
      info = name
    )
  })
})

test_that("screen_report prints a readable report", {
  screen <- make_screen_for_report(
    class = "MultiplexScreen",
    failed_queries = paste0("Q", seq_len(10L))
  )

  output <- capture.output(screen_report(screen, print = TRUE, width = 12))
  text <- paste(output, collapse = "\n")
  markers <- c(
    "CeRberus screen report",
    "OVERVIEW",
    "MODEL",
    "CHECKS",
    "PROBLEMS",
    "RESULTS",
    "query_sufficient: OK",
    "library_sufficient: PROBLEM",
    "+2 more"
  )

  expect_equal(output[[2L]], paste(rep("=", 12L), collapse = ""))
  expect_true(all(vapply(markers, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("screen_report writes YAML variants without printing", {
  screen <- make_screen_for_report(class = "MultiplexScreen")
  directory <- tempfile("screen-report-directory-")

  for (extension in c(".yaml", ".yml")) {
    path <- file.path(directory, "nested", paste0("screen-report", extension))

    output <- capture.output(
      report <- screen_report(screen, file = path)
    )
    written <- yaml::read_yaml(path)

    expect_equal(
      list(
        output = output,
        exists = file.exists(path),
        report = report[c("report_version", "screen", "dimensions", "checks")],
        written = written[c("report_version", "screen", "dimensions", "checks")]
      ),
      list(
        output = character(),
        exists = TRUE,
        report = report[c("report_version", "screen", "dimensions", "checks")],
        written = report[c("report_version", "screen", "dimensions", "checks")]
      ),
      info = extension
    )
  }
})

test_that("screen_report validates output controls", {
  screen <- make_screen_for_report(class = "MultiplexScreen")
  cases <- list(
    file = list(
      values = list(NA_character_, character(), c("a.yaml", "b.yaml"), 1),
      call = function(value) screen_report(screen, file = value, print = FALSE),
      error = "file must be NULL or a single path string"
    ),
    extension = list(
      values = list("report.txt"),
      call = function(value) screen_report(screen, file = value, print = FALSE),
      error = "file must use a .yaml or .yml extension"
    ),
    print = list(
      values = list(NA, logical(), c(TRUE, FALSE), 1, "TRUE"),
      call = function(value) screen_report(screen, print = value),
      error = "print must be TRUE or FALSE"
    ),
    width = list(
      values = list(NA_real_, Inf, 0, -1, numeric(), c(40, 80), "80"),
      call = function(value) {
        screen_report(screen, print = FALSE, width = value)
      },
      error = "width must be one finite positive numeric value"
    )
  )

  purrr::iwalk(cases, function(case, name) {
    purrr::walk(case$values, function(value) {
      expect_error(case$call(value), case$error, info = name)
    })
  })
})
