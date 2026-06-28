make_configuration_screen <- function(dupcor) {
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
    metadata = list(),
    checks = list(),
    errors = list()
  )
}

test_that("find_optimal_configuration selects the weakest positive usable correlation", {
  screens <- list(
    negative = make_configuration_screen(-0.25),
    zero = make_configuration_screen(0),
    low_positive = make_configuration_screen(c(0.041, 0.049)),
    high_positive = make_configuration_screen(0.3),
    missing = make_configuration_screen(c(NA_real_, NaN))
  )

  result <- CeRberus:::find_optimal_configuration(screens)

  expect_type(result, "list")
  expect_named(result, "low_positive")
  expect_s4_class(result[[1L]], "ScreenBase")

  dupcor_data <- result[[1L]]@metadata$dupcor_data
  expect_s3_class(dupcor_data, "data.table")
  expect_equal(
    dupcor_data$config,
    c("negative", "zero", "low_positive", "high_positive", "missing")
  )
  expect_equal(dupcor_data$dcor, c(-0.25, 0, 0.045, 0.3, NA))
  expect_equal(dupcor_data$kept, c("", "", "selected", "", ""))
})

test_that("find_optimal_configuration returns all annotated configurations when keep_all is TRUE", {
  screens <- list(
    first = make_configuration_screen(0.2),
    selected = make_configuration_screen(0.1),
    third = make_configuration_screen(0.4)
  )

  result <- CeRberus:::find_optimal_configuration(screens, keep_all = TRUE)

  expect_named(result, names(screens))
  expect_equal(length(result), 3L)

  for (screen in result) {
    expect_s3_class(screen@metadata$dupcor_data, "data.table")
    expect_equal(screen@metadata$dupcor_data$config, names(screens))
    expect_equal(screen@metadata$dupcor_data$kept, c("", "selected", ""))
  }
})

test_that("find_optimal_configuration selects zero when no positive correlation exists", {
  screens <- list(
    negative = make_configuration_screen(-0.4),
    zero = make_configuration_screen(0),
    more_negative = make_configuration_screen(-0.8)
  )

  result <- CeRberus:::find_optimal_configuration(screens)

  expect_named(result, "zero")
  expect_equal(result[[1L]]@metadata$dupcor_data$kept, c("", "selected", ""))
})

test_that("find_optimal_configuration selects the least negative correlation when all usable correlations are negative", {
  screens <- list(
    very_negative = make_configuration_screen(-0.7),
    least_negative = make_configuration_screen(-0.1),
    moderately_negative = make_configuration_screen(-0.3)
  )

  result <- CeRberus:::find_optimal_configuration(screens)

  expect_named(result, "least_negative")
  expect_equal(result[[1L]]@metadata$dupcor_data$kept, c("", "selected", ""))
})

test_that("find_optimal_configuration averages duplicate-correlation vectors before selection", {
  screens <- list(
    selected = make_configuration_screen(c(0.05, 0.15, NA_real_)),
    not_selected = make_configuration_screen(c(0.2, 0.4))
  )

  result <- CeRberus:::find_optimal_configuration(screens)

  expect_named(result, "selected")
  expect_equal(result[[1L]]@metadata$dupcor_data$dcor, c(0.1, 0.3))
})

test_that("find_optimal_configuration validates input list and keep_all", {
  valid_screen <- make_configuration_screen(0.1)

  expect_error(
    CeRberus:::find_optimal_configuration(list()),
    "GI_list must contain at least one screen object"
  )
  expect_error(
    CeRberus:::find_optimal_configuration(list(valid_screen)),
    "GI_list must be named"
  )
  expect_error(
    CeRberus:::find_optimal_configuration(setNames(list(valid_screen), "")),
    "GI_list must be named"
  )
  expect_error(
    CeRberus:::find_optimal_configuration(list(a = valid_screen, a = valid_screen)),
    "GI_list must have unique names"
  )
  expect_error(
    CeRberus:::find_optimal_configuration(list(a = valid_screen), keep_all = NA),
    "keep_all must be TRUE or FALSE"
  )
  expect_error(
    CeRberus:::find_optimal_configuration(list(a = valid_screen), keep_all = c(TRUE, FALSE)),
    "keep_all must be TRUE or FALSE"
  )
})

test_that("find_optimal_configuration errors when no finite duplicate-correlation estimates are available", {
  screens <- list(
    missing = make_configuration_screen(c(NA_real_, NaN)),
    infinite = make_configuration_screen(Inf)
  )

  expect_error(
    CeRberus:::find_optimal_configuration(screens),
    "No usable duplicate-correlation estimates found"
  )
})
