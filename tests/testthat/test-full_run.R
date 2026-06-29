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
    screen_attr = list(),
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
  calls = new.env(parent = emptyenv())
) {
  calls$collected_input <- NULL
  calls$make_pos_agnostic <- NULL
  calls$collect_verbose <- NULL
  calls$plot_path <- NULL
  calls$plot_verbose <- NULL
  calls$fdr_method <- NULL
  calls$keep_all <- NULL
  calls$screen_report_called <- FALSE

  testthat::local_mocked_bindings(
    collect_all_layer_configurations = function(
      GI_data,
      make_pos_agnostic,
      verbose = FALSE
    ) {
      calls$collected_input <- as.data.frame(GI_data)
      calls$make_pos_agnostic <- make_pos_agnostic
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
    compute_dupCorrelation = function(.x, ...) {
      .x@metadata$compute_dupCorrelation_called <- TRUE
      .x
    },
    find_optimal_configuration = function(GI_list, keep_all = FALSE) {
      calls$keep_all <- keep_all
      if (isTRUE(keep_all)) {
        return(GI_list)
      }
      GI_list["default_guide_pair_used"]
    },
    compute_dupcor_plot = function(.data, .fpath, verbose = FALSE) {
      calls$plot_path <- .fpath
      calls$plot_verbose <- verbose
      dir.create(dirname(.fpath), showWarnings = FALSE, recursive = TRUE)
      file.create(.fpath)
      .data
    },
    compute_models = function(GI_obj, ...) {
      GI_obj@metadata$compute_models_called <- TRUE
      GI_obj
    },
    collect_GIs = function(GI_obj, FDR_method, ...) {
      calls$fdr_method <- FDR_method
      GI_obj@metadata$collect_GIs_called <- TRUE
      GI_obj
    },
    GI_df = function(GI_obj, ...) {
      data.table::data.table(
        gene_pair = GI_obj@metadata$config,
        GI = GI_obj@dupCorrelation
      )
    },
    screenReport = function(GI_obj, interactive = FALSE, print = TRUE, ...) {
      calls$screen_report_called <- TRUE
      invisible(list())
    },
    .package = "CeRberus"
  )

  force(code)
}

test_that("full_run reads CSV scores and forwards options to the pipeline", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  data.table::fwrite(make_full_run_scores(), scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    FDR = "bonferroni",
    overwrite_output = FALSE,
    verbose = TRUE,
    make_symmetric = TRUE
  )
  calls <- new.env(parent = emptyenv())

  result <- with_mocked_full_run_pipeline(
    full_run(yaml_fpath),
    calls = calls
  )

  expect_type(result, "list")
  expect_named(result, "default_guide_pair_used")
  expect_s4_class(result[[1]], "ScreenBase")
  expect_equal(calls$collected_input, make_full_run_scores())
  expect_true(calls$make_pos_agnostic)
  expect_true(calls$collect_verbose)
  expect_equal(calls$fdr_method, "bonferroni")
  expect_null(calls$plot_verbose)
  expect_null(calls$plot_path)
  expect_true(calls$screen_report_called)
})

test_that("full_run reads RDS score files", {
  scores_file <- tempfile(fileext = ".rds")
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  scores <- make_full_run_scores()
  saveRDS(scores, scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = FALSE
  )
  calls <- new.env(parent = emptyenv())

  result <- with_mocked_full_run_pipeline(
    full_run(yaml_fpath),
    calls = calls
  )

  expect_named(result, "default_guide_pair_used")
  expect_equal(calls$collected_input, scores)
})

test_that("full_run creates the configured output directory", {
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

  with_mocked_full_run_pipeline(
    full_run(yaml_fpath)
  )

  expect_true(dir.exists(output_directory))
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

  expect_true(file.exists(file.path(output_directory, "all_GI_objects.rds")))
  expect_true(file.exists(file.path(
    output_directory,
    "duplicateCorrelationPlot.png"
  )))
  expect_true(file.exists(file.path(
    output_directory,
    "duplicate_correlation.csv"
  )))
  expect_true(file.exists(file.path(
    output_directory,
    "GI_scores_default_guide_pair_used.csv"
  )))
  expect_equal(
    calls$plot_path,
    file.path(
      normalizePath(output_directory, winslash = "/", mustWork = FALSE),
      "duplicateCorrelationPlot.png"
    )
  )

  intermediate <- readRDS(file.path(output_directory, "all_GI_objects.rds"))
  expect_named(
    intermediate,
    c("default_guide_pair_used", "default_tech_rep_used")
  )
  expect_named(result, "default_guide_pair_used")
})

test_that("full_run does not write CSV/RDS outputs when overwrite_output is FALSE", {
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

  with_mocked_full_run_pipeline(
    full_run(yaml_fpath)
  )

  expect_false(file.exists(file.path(output_directory, "all_GI_objects.rds")))
  expect_false(file.exists(file.path(
    output_directory,
    "duplicate_correlation.csv"
  )))
  expect_false(file.exists(file.path(
    output_directory,
    "GI_scores_default_guide_pair_used.csv"
  )))
  expect_false(file.exists(file.path(
    output_directory,
    "duplicateCorrelationPlot.png"
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
    overwrite_output = FALSE,
    keep_all_configurations = TRUE
  )
  calls <- new.env(parent = emptyenv())

  result <- with_mocked_full_run_pipeline(
    full_run(yaml_fpath),
    calls = calls
  )

  expect_true(calls$keep_all)
  expect_named(result, c("default_guide_pair_used", "default_tech_rep_used"))
})

test_that("full_run rejects scores files without an extension", {
  scores_file <- tempfile()
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  file.create(scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory
  )

  expect_error(
    full_run(yaml_fpath),
    "scores_file must have a file extension"
  )
})

test_that("full_run rejects unsupported scores file extensions", {
  scores_file <- tempfile(fileext = ".txt")
  output_directory <- tempfile("full-run-output-")
  yaml_fpath <- tempfile(fileext = ".yaml")
  writeLines("not,a,supported,file", scores_file)
  write_full_run_instructions(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory
  )

  expect_error(
    full_run(yaml_fpath),
    "Unsupported scores_file extension: \\.txt"
  )
})
