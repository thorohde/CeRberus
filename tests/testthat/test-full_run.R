write_full_run_instructions <- function(
  path,
  scores_file,
  output_directory,
  ...
) {
  yaml::write_yaml(
    c(
      list(
        scores_file = scores_file,
        output_directory = output_directory
      ),
      list(...)
    ),
    file = path
  )
}

make_full_run_scores <- function() {
  data.frame(
    query_gene = c("A", "A", "B", "B"),
    library_gene = c("C", "D", "C", "D"),
    bio_rep = c("b1", "b1", "b2", "b2"),
    tech_rep = c("t1", "t2", "t1", "t2"),
    guide_pair = c("g1", "g2", "g1", "g2"),
    GI = c(0.1, 0.2, -0.1, -0.2)
  )
}

make_full_run_screen <- function(
  name = "default_guide_pair_used",
  dupcor = 0.1
) {
  methods::new(
    "ScreenBase",
    guideLFCs = methods::new(
      "gRNA_LFC",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character()
    ),
    guideGIs = methods::new(
      "gRNA_GI",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character(),
      block_layer = character(),
      blocks = character(),
      use_blocks = FALSE,
      block_description = character(),
      collapse = character()
    ),
    limma_models = list(),
    geneGIs = array(numeric(), dim = 0),
    screen_attr = methods::new("ScreenDesign"),
    dupCorrelation = dupcor,
    metadata = list(
      config = name,
      dupcor_data = data.table::data.table(
        config = name,
        dcor = dupcor,
        kept = "selected"
      )
    ),
    checks = list(),
    errors = list()
  )
}

with_mocked_full_run_pipeline <- function(
  code,
  calls = new.env(parent = emptyenv()),
  inspect_output_directory = NULL
) {
  calls$collected_input <- NULL
  calls$screen_type <- NULL
  calls$pos_agnostic <- NULL
  calls$symmetric_analysis_method <- NULL
  calls$non_targeting_controls <- NULL
  calls$retain_directional <- NULL
  calls$collect_verbose <- NULL
  calls$plot_path <- NULL
  calls$plot_verbose <- NULL
  calls$qq_plot_paths <- character()
  calls$fdr_method <- NULL
  calls$keep_all <- NULL
  calls$screen_report_called <- FALSE
  calls$output_contents_at_collection <- NULL

  testthat::local_mocked_bindings(
    collect_all_layer_configurations = function(
      gi_data,
      screen_type = "auto",
      pos_agnostic,
      symmetric_analysis_method = "preaverage",
      non_targeting_controls = NULL,
      retain_directional = FALSE,
      verbose = FALSE
    ) {
      if (!is.null(inspect_output_directory)) {
        calls$output_contents_at_collection <- list.files(
          inspect_output_directory,
          all.files = TRUE,
          no.. = TRUE
        )
      }
      calls$collected_input <- as.data.frame(gi_data)
      calls$screen_type <- screen_type
      calls$pos_agnostic <- pos_agnostic
      calls$symmetric_analysis_method <- symmetric_analysis_method
      calls$non_targeting_controls <- non_targeting_controls
      calls$retain_directional <- retain_directional
      calls$collect_verbose <- verbose
      list(
        default_guide_pair_used = make_full_run_screen(
          "default_guide_pair_used",
          0.2
        ),
        default_tech_rep_used = make_full_run_screen(
          "default_tech_rep_used",
          0.4
        )
      )
    },
    compute_dup_correlation = function(.x, ...) {
      .x@metadata$compute_dup_correlation_called <- TRUE
      .x
    },
    find_optimal_configuration = function(gi_list, keep_all = FALSE) {
      calls$keep_all <- keep_all
      dupcor_data <- data.table::data.table(
        config = names(gi_list),
        dcor = c(0.2, 0.4),
        kept = c("selected", "")
      )
      for (name in names(gi_list)) {
        gi_list[[name]]@metadata$dupcor_data <- dupcor_data
      }
      if (isTRUE(keep_all)) {
        return(gi_list)
      }
      gi_list["default_guide_pair_used"]
    },
    compute_dupcor_plot = function(.data, .fpath, verbose = FALSE) {
      calls$plot_path <- .fpath
      calls$plot_verbose <- verbose
      dir.create(dirname(.fpath), showWarnings = FALSE, recursive = TRUE)
      file.create(.fpath)
      .data
    },
    pvalue_qq_plot = function(gi_obj, .fpath, verbose = FALSE) {
      calls$qq_plot_paths <- c(calls$qq_plot_paths, .fpath)
      dir.create(dirname(.fpath), showWarnings = FALSE, recursive = TRUE)
      file.create(.fpath)
      gi_obj
    },
    compute_models = function(gi_obj, ...) {
      gi_obj@metadata$compute_models_called <- TRUE
      gi_obj
    },
    collect_gis = function(gi_obj, fdr_method, ...) {
      calls$fdr_method <- fdr_method
      gi_obj@metadata$collect_gis_called <- TRUE
      gi_obj
    },
    gi_df = function(gi_obj, ...) {
      data.table::data.table(
        gene_pair = gi_obj@metadata$config,
        GI = gi_obj@dupCorrelation
      )
    },
    screen_report = function(
      gi_obj,
      file = NULL,
      print = is.null(file),
      ...
    ) {
      calls$screen_report_called <- TRUE
      report <- list(
        report_version = "1.1",
        screen = list(configuration = gi_obj@metadata$config)
      )

      report
    },
    .package = "CeRberus"
  )

  force(code)
}

test_that("full_run reads supported score files and forwards pipeline options", {
  writers <- list(
    csv = function(scores, path) data.table::fwrite(scores, path),
    rds = function(scores, path) saveRDS(scores, path)
  )

  purrr::iwalk(writers, function(write_scores, extension) {
    scores_file <- tempfile(fileext = paste0(".", extension))
    output_directory <- tempfile("full-run-output-")
    yaml_fpath <- tempfile(fileext = ".yaml")
    scores <- make_full_run_scores()
    write_scores(scores, scores_file)
    write_full_run_instructions(
      yaml_fpath,
      scores_file = scores_file,
      output_directory = output_directory,
      FDR = "bonferroni",
      overwrite_output = FALSE,
      verbose = TRUE,
      screen_type = "fixed_pair",
      pos_agnostic = TRUE,
      symmetric_analysis_method = "preaverage",
      retain_directional = TRUE
    )
    calls <- new.env(parent = emptyenv())

    result <- with_mocked_full_run_pipeline(
      full_run(yaml_fpath),
      calls = calls
    )

    expect_type(result, "list")
    expect_named(result, "default_guide_pair_used", info = extension)
    expect_s4_class(result[[1]], "ScreenBase")
    expect_equal(calls$collected_input, scores, info = extension)
    expect_identical(calls$screen_type, "fixed_pair", info = extension)
    expect_true(calls$pos_agnostic, info = extension)
    expect_identical(
      calls$symmetric_analysis_method,
      "preaverage",
      info = extension
    )
    expect_true(calls$retain_directional, info = extension)
    expect_true(calls$collect_verbose, info = extension)
    expect_equal(calls$fdr_method, "bonferroni", info = extension)
    expect_null(calls$plot_verbose, info = extension)
    expect_null(calls$plot_path, info = extension)
    expect_true(calls$screen_report_called, info = extension)
  })
})

test_that("full_run forwards YAML non-targeting controls", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  data.table::fwrite(make_full_run_scores(), scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = FALSE,
    non_targeting_controls = c("NTC_1", "NTC_2")
  )
  calls <- new.env(parent = emptyenv())

  with_mocked_full_run_pipeline(full_run(yaml_fpath), calls = calls)

  expect_identical(calls$non_targeting_controls, c("NTC_1", "NTC_2"))
})

test_that("full_run preserves output directories when overwrite is disabled", {
  artifacts <- c(
    "all_gi_objects.rds",
    "duplicate_correlation.csv",
    "GI_scores_default_guide_pair_used.csv",
    "duplicateCorrelationPlot.png",
    "pvalueQQPlot_default_guide_pair_used.png",
    "screen_report.yaml",
    "CeRberus.log"
  )

  purrr::walk(c(FALSE, TRUE), function(directory_exists) {
    scores_file <- tempfile(fileext = ".csv")
    output_directory <- tempfile("full-run-output-")
    yaml_fpath <- tempfile(fileext = ".yaml")
    stale_file <- file.path(output_directory, "stale.csv")
    if (directory_exists) {
      dir.create(output_directory)
      file.create(stale_file)
    }
    data.table::fwrite(make_full_run_scores(), scores_file)
    write_full_run_instructions(
      yaml_fpath,
      scores_file = scores_file,
      output_directory = output_directory,
      overwrite_output = FALSE
    )

    with_mocked_full_run_pipeline(full_run(yaml_fpath))

    case_name <- if (directory_exists) "existing" else "missing"
    expect_true(dir.exists(output_directory), info = case_name)
    purrr::walk(artifacts, function(artifact) {
      expect_false(
        file.exists(file.path(output_directory, artifact)),
        info = paste(case_name, artifact)
      )
    })
    if (directory_exists) {
      expect_true(file.exists(stale_file), info = case_name)
    }
  })
})

test_that("full_run removes previous outputs and preserves other output files", {
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- file.path(output_directory, "instructions.yaml")
  scores_file <- file.path(output_directory, "scores.csv")
  dir.create(file.path(output_directory, "nested"), recursive = TRUE)
  file.create(
    file.path(output_directory, "stale.csv"),
    file.path(output_directory, ".hidden"),
    file.path(output_directory, "nested", "stale.txt"),
    file.path(output_directory, "CeRberus_error.log")
  )
  data.table::fwrite(make_full_run_scores(), scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = TRUE
  )
  calls <- new.env(parent = emptyenv())

  with_mocked_full_run_pipeline(
    full_run(yaml_fpath),
    calls = calls,
    inspect_output_directory = output_directory
  )

  expect_setequal(
    calls$output_contents_at_collection,
    c("scores.csv", "instructions.yaml", "nested", "stale.csv", ".hidden")
  )
  expect_true(file.exists(scores_file))
  expect_true(file.exists(yaml_fpath))
  expect_true(file.exists(file.path(output_directory, "stale.csv")))
  expect_true(file.exists(file.path(output_directory, ".hidden")))
  expect_true(file.exists(file.path(output_directory, "nested", "stale.txt")))
  expect_false(file.exists(file.path(output_directory, "CeRberus_error.log")))
  expect_true(file.exists(file.path(
    output_directory,
    "GI_scores_default_guide_pair_used.csv"
  )))
})

test_that("full_run refuses to empty a filesystem root", {
  scores_file <- tempfile(fileext = ".csv")
  yaml_fpath <- tempfile(fileext = ".yaml")
  root_directory <- normalizePath(
    .Platform$file.sep,
    winslash = "/",
    mustWork = TRUE
  )
  data.table::fwrite(make_full_run_scores(), scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = root_directory,
    overwrite_output = TRUE
  )

  suppressWarnings(
    expect_error(
      full_run(yaml_fpath),
      "Refusing to empty a filesystem root",
      fixed = TRUE
    )
  )
})

test_that("full_run rejects an output path that is a regular file", {
  scores_file <- tempfile(fileext = ".csv")
  output_path <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  data.table::fwrite(make_full_run_scores(), scores_file)
  file.create(output_path)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_path,
    overwrite_output = TRUE
  )

  expect_error(
    full_run(yaml_fpath),
    "output_directory exists but is not a directory",
    fixed = TRUE
  )
  expect_true(file.exists(output_path))
  expect_false(dir.exists(output_path))
})

test_that("full_run writes intermediate and final outputs when overwrite_output is TRUE", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  data.table::fwrite(make_full_run_scores(), scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = TRUE
  )
  calls <- new.env(parent = emptyenv())

  result <- with_mocked_full_run_pipeline(
    full_run(yaml_fpath),
    calls = calls
  )

  expect_true(file.exists(file.path(output_directory, "all_gi_objects.rds")))
  expect_true(file.exists(file.path(
    output_directory,
    "duplicateCorrelationPlot.png"
  )))
  expect_true(file.exists(file.path(
    output_directory,
    "pvalueQQPlot_default_guide_pair_used.png"
  )))
  expect_true(file.exists(file.path(
    output_directory,
    "duplicate_correlation.csv"
  )))
  expect_true(file.exists(file.path(
    output_directory,
    "GI_scores_default_guide_pair_used.csv"
  )))
  report_path <- file.path(
    normalizePath(output_directory, winslash = "/", mustWork = FALSE),
    "screen_report.yaml"
  )
  expect_true(file.exists(report_path))
  report <- yaml::read_yaml(report_path)
  expect_identical(report$report_version, "1.1")
  expect_identical(
    report$selection$selected_configuration,
    "default_guide_pair_used"
  )
  expect_named(
    report$selection$evaluated_configurations,
    c("default_guide_pair_used", "default_tech_rep_used")
  )
  expect_true(
    report$selection$evaluated_configurations$default_guide_pair_used$selected
  )
  expect_false(
    report$selection$evaluated_configurations$default_tech_rep_used$selected
  )
  expect_named(report$configurations, "default_guide_pair_used")
  expect_identical(
    report$configurations$default_guide_pair_used$screen$configuration,
    "default_guide_pair_used"
  )
  expect_true(file.exists(file.path(output_directory, "CeRberus.log")))
  success_log <- paste(
    readLines(file.path(output_directory, "CeRberus.log"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(success_log, "Status: completed", fixed = TRUE)
  expect_match(success_log, "Pipeline stage: complete", fixed = TRUE)
  expect_match(
    success_log,
    "Configuration entry: default_guide_pair_used",
    fixed = TRUE
  )
  expect_equal(
    calls$plot_path,
    file.path(
      normalizePath(output_directory, winslash = "/", mustWork = FALSE),
      "duplicateCorrelationPlot.png"
    )
  )

  intermediate <- readRDS(file.path(output_directory, "all_gi_objects.rds"))
  expect_named(
    intermediate,
    c("default_guide_pair_used", "default_tech_rep_used")
  )
  expect_named(result, "default_guide_pair_used")
  expect_equal(
    calls$qq_plot_paths,
    file.path(
      normalizePath(output_directory, winslash = "/", mustWork = FALSE),
      "pvalueQQPlot_default_guide_pair_used.png"
    )
  )
})

test_that("full_run writes a failure log and preserves the original error", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  data.table::fwrite(make_full_run_scores(), scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = FALSE
  )

  testthat::local_mocked_bindings(
    collect_all_layer_configurations = function(...) {
      list(
        default_guide_pair_used = make_full_run_screen(
          "default_guide_pair_used",
          0.2
        )
      )
    },
    compute_dup_correlation = function(.x, ...) .x,
    find_optimal_configuration = function(gi_list, ...) gi_list,
    compute_models = function(gi_obj, ...) {
      stop("deliberate model failure", call. = FALSE)
    },
    .package = "CeRberus"
  )

  expect_error(
    full_run(yaml_fpath),
    "deliberate model failure",
    fixed = TRUE
  )

  error_log_path <- file.path(output_directory, "CeRberus_error.log")
  expect_true(file.exists(error_log_path))
  error_log <- paste(readLines(error_log_path, warn = FALSE), collapse = "\n")
  expect_match(error_log, "Status: failed", fixed = TRUE)
  expect_match(error_log, "Pipeline stage: fit_models", fixed = TRUE)
  expect_match(error_log, "Condition: i In index: 1.", fixed = TRUE)
  expect_match(error_log, "deliberate model failure", fixed = TRUE)
  expect_match(
    error_log,
    "Configuration entry: default_guide_pair_used",
    fixed = TRUE
  )
  expect_false(file.exists(file.path(output_directory, "CeRberus.log")))
  expect_false(file.exists(file.path(
    output_directory,
    "screen_report.yaml"
  )))
})

test_that("full_run can return NULL after running the pipeline", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  data.table::fwrite(make_full_run_scores(), scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = TRUE
  )

  result <- with_mocked_full_run_pipeline(
    full_run(yaml_fpath, return_output = FALSE)
  )

  expect_null(result)
  expect_true(file.exists(file.path(
    output_directory,
    "GI_scores_default_guide_pair_used.csv"
  )))
  expect_true(file.exists(file.path(
    output_directory,
    "screen_report.yaml"
  )))
})

test_that("full_run forwards keep_all_configurations to configuration selection", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  data.table::fwrite(make_full_run_scores(), scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = TRUE,
    keep_all_configurations = TRUE
  )
  calls <- new.env(parent = emptyenv())

  result <- with_mocked_full_run_pipeline(
    full_run(yaml_fpath),
    calls = calls
  )

  expect_true(calls$keep_all)
  expect_named(result, c("default_guide_pair_used", "default_tech_rep_used"))
  expect_equal(
    calls$qq_plot_paths,
    file.path(
      normalizePath(output_directory, winslash = "/", mustWork = FALSE),
      paste0(
        "pvalueQQPlot_",
        c("default_guide_pair_used", "default_tech_rep_used"),
        ".png"
      )
    )
  )
  report_path <- file.path(
    output_directory,
    "screen_report.yaml"
  )
  expect_true(file.exists(report_path))
  report <- yaml::read_yaml(report_path)
  expect_named(
    report$configurations,
    c("default_guide_pair_used", "default_tech_rep_used")
  )
  expect_identical(
    report$configurations$default_tech_rep_used$screen$configuration,
    "default_tech_rep_used"
  )
})

test_that("full_run rejects invalid scores file extensions", {
  cases <- list(
    missing = list(
      fileext = "",
      error = "scores_file must have a file extension"
    ),
    unsupported = list(
      fileext = ".txt",
      error = "Unsupported scores_file extension: \\.txt"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    scores_file <- tempfile(fileext = case$fileext)
    output_directory <- tempfile("full-run-output-")
    yaml_fpath <- tempfile(fileext = ".yaml")
    file.create(scores_file)
    write_full_run_instructions(
      yaml_fpath,
      scores_file = scores_file,
      output_directory = output_directory
    )

    expect_error(full_run(yaml_fpath), case$error, info = case_name)
  })
})
