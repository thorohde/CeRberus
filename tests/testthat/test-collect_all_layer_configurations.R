make_layer_configuration_scores <- function(include = c("tech_rep", "bio_rep", "guide_pair")) {
  base <- data.frame(
    query_gene = c("A", "A", "B", "B"),
    library_gene = c("A", "B", "A", "B"),
    GI = c(0.1, 0.2, -0.1, -0.2)
  )

  if ("tech_rep" %in% include) {
    base$tech_rep <- c("t1", "t2", "t1", "t2")
  }
  if ("bio_rep" %in% include) {
    base$bio_rep <- c("b1", "b1", "b2", "b2")
  }
  if ("guide_pair" %in% include) {
    base$guide_pair <- c("g1", "g2", "g1", "g2")
  }

  base
}

with_mocked_GIScores <- function(code, calls = new.env(parent = emptyenv())) {
  calls$args <- list()

  testthat::local_mocked_bindings(
    GIScores = function(input,
                        collapse_layers = NULL,
                        block_layer = NULL,
                        pos_agnostic = FALSE,
                        verbose = FALSE,
                        ...) {
      calls$args[[length(calls$args) + 1L]] <- list(
        input = input,
        collapse_layers = collapse_layers,
        block_layer = block_layer,
        pos_agnostic = pos_agnostic,
        verbose = verbose
      )

      list(
        collapse_layers = collapse_layers,
        block_layer = block_layer,
        pos_agnostic = pos_agnostic,
        verbose = verbose
      )
    },
    .package = "CeRberus"
  )

  force(code)
}

test_that("collect_all_layer_configurations builds all default and collapsed layer configurations", {
  scores <- make_layer_configuration_scores()
  calls <- new.env(parent = emptyenv())

  result <- with_mocked_GIScores(
    CeRberus:::collect_all_layer_configurations(
      scores,
      make_pos_agnostic = TRUE,
      verbose = TRUE
    ),
    calls = calls
  )

  expect_type(result, "list")
  expect_length(result, 12L)
  expect_named(
    result,
    c(
      "default_tech_rep_used",
      "default_bio_rep_used",
      "default_guide_pair_used",
      "tech_rep_collapsed_bio_rep_used",
      "tech_rep_collapsed_guide_pair_used",
      "bio_rep_collapsed_tech_rep_used",
      "bio_rep_collapsed_guide_pair_used",
      "guide_pair_collapsed_tech_rep_used",
      "guide_pair_collapsed_bio_rep_used",
      "tech_rep_bio_rep_collapsed_guide_pair_used",
      "tech_rep_guide_pair_collapsed_bio_rep_used",
      "bio_rep_guide_pair_collapsed_tech_rep_used"
    )
  )

  expect_length(calls$args, 12L)
  expect_true(all(purrr::map_lgl(calls$args, "pos_agnostic")))
  expect_true(all(purrr::map_lgl(calls$args, "verbose")))
  expect_true(all(purrr::map_lgl(calls$args, ~ identical(.x$input, scores))))

  expect_null(calls$args[[1L]]$collapse_layers)
  expect_equal(calls$args[[1L]]$block_layer, "tech_rep")
  expect_equal(calls$args[[4L]]$collapse_layers, "tech_rep")
  expect_equal(calls$args[[4L]]$block_layer, "bio_rep")
  expect_equal(calls$args[[10L]]$collapse_layers, c("tech_rep", "bio_rep"))
  expect_equal(calls$args[[10L]]$block_layer, "guide_pair")
})

test_that("collect_all_layer_configurations respects requested use layers and ignores absent columns", {
  scores <- make_layer_configuration_scores(include = c("tech_rep", "bio_rep"))
  calls <- new.env(parent = emptyenv())

  result <- with_mocked_GIScores(
    CeRberus:::collect_all_layer_configurations(
      scores,
      .to_use = c("bio_rep", "guide_pair", "not_a_column"),
      make_pos_agnostic = FALSE,
      verbose = FALSE
    ),
    calls = calls
  )

  expect_named(
    result,
    c(
      "default_bio_rep_used",
      "tech_rep_collapsed_bio_rep_used",
      "bio_rep_collapsed_tech_rep_used"
    )
  )
  expect_length(calls$args, 3L)
  expect_equal(purrr::map_chr(calls$args, "block_layer"), c("bio_rep", "bio_rep", "tech_rep"))
  expect_equal(calls$args[[1L]]$collapse_layers, NULL)
  expect_equal(calls$args[[2L]]$collapse_layers, "tech_rep")
  expect_equal(calls$args[[3L]]$collapse_layers, "bio_rep")
  expect_false(any(purrr::map_lgl(calls$args, "pos_agnostic")))
  expect_false(any(purrr::map_lgl(calls$args, "verbose")))
})

test_that("collect_all_layer_configurations only creates defaults when a single collapsible layer is present", {
  scores <- make_layer_configuration_scores(include = "guide_pair")
  calls <- new.env(parent = emptyenv())

  result <- with_mocked_GIScores(
    CeRberus:::collect_all_layer_configurations(
      scores,
      make_pos_agnostic = FALSE
    ),
    calls = calls
  )

  expect_named(result, "default_guide_pair_used")
  expect_length(calls$args, 1L)
  expect_equal(calls$args[[1L]]$block_layer, "guide_pair")
  expect_null(calls$args[[1L]]$collapse_layers)
})