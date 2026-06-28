make_screen_for_dupCorrelation_df <- function(dupcor) {
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

test_that("dupCorrelation_df returns a one-column data.table for scalar duplicate correlation", {
  screen <- make_screen_for_dupCorrelation_df(0.123)

  result <- dupCorrelation_df(screen)

  expect_s3_class(result, "data.table")
  expect_named(result, "dupcor")
  expect_equal(nrow(result), 1L)
  expect_equal(result$dupcor, 0.123)
})

test_that("dupCorrelation_df preserves vector duplicate correlations", {
  screen <- make_screen_for_dupCorrelation_df(c(Q1 = 0.1, Q2 = 0.2, Q3 = NA_real_))

  result <- dupCorrelation_df(screen)

  expect_s3_class(result, "data.table")
  expect_named(result, "dupcor")
  expect_equal(nrow(result), 3L)
  expect_equal(result$dupcor, c(0.1, 0.2, NA_real_))
})

test_that("dupCorrelation_df handles empty duplicate-correlation vectors", {
  screen <- make_screen_for_dupCorrelation_df(numeric())

  result <- dupCorrelation_df(screen)

  expect_s3_class(result, "data.table")
  expect_named(result, "dupcor")
  expect_equal(nrow(result), 0L)
  expect_type(result$dupcor, "double")
})

test_that("dupCorrelation_df does not modify the input screen", {
  screen <- make_screen_for_dupCorrelation_df(c(0.1, 0.2))
  original <- screen

  dupCorrelation_df(screen)

  expect_equal(screen, original)
})

test_that("dupcor_df alias generic has no ScreenBase method implemented", {
  screen <- make_screen_for_dupCorrelation_df(0.1)

  expect_error(
    dupcor_df(screen),
    "unable to find an inherited method"
  )
})
