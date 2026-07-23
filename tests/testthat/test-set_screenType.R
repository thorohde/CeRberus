make_multiplex_like_scores_for_set_screen_type <- function() {
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

make_fixed_pair_like_scores_for_set_screen_type <- function() {
  data.frame(
    query_gene = rep(c("A", "B"), each = 4),
    library_gene = rep(c("C", "D"), each = 4),
    bio_rep = rep(c("b1", "b2"), 4),
    tech_rep = rep(c("t1", "t2"), each = 2, times = 2),
    guide_pair = rep(c("g1", "g2"), 4),
    GI = seq_len(8)
  )
}

test_that("set_screen_type keeps multiplex-compatible screens as multiplex", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = TRUE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE
  )

  result <- set_screen_type(as(screen, "ScreenBase"))

  expect_s4_class(result, "MultiplexScreen")
  expect_equal(result@guideGIs@space, c("query_gene", "library_gene"))
  expect_equal(result@guideGIs@replicates, "guide_pair")
  expect_equal(result@metadata$input, screen@metadata$input)
})

test_that("set_screen_type switches to fixed-pair when library checks fail", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = FALSE,
    stable_library_size = FALSE,
    sufficient_tests_per_query = FALSE
  )

  result <- set_screen_type(as(screen, "ScreenBase"))

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
  expect_equal(result@guideGIs@replicates, "guide_pair")
})

test_that("set_screen_type keeps fixed-pair fallback when multiplex criteria are not met", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = FALSE,
    library_sufficient = TRUE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE
  )

  result <- set_screen_type(as(screen, "ScreenBase"))

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
})

test_that("set_screen_type honors requested fixed-pair type when checks allow multiplex", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = TRUE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE
  )
  screen@metadata$requested_screen_type <- "fixed_pair"

  expect_warning(
    result <- set_screen_type(as(screen, "ScreenBase")),
    "overrides inferred screen type 'multiplex'"
  )

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
  expect_equal(result@guideGIs@replicates, "guide_pair")
  expect_identical(result@metadata$inferred_screen_type, "multiplex")
  expect_identical(result@metadata$selected_screen_type, "fixed_pair")
})

test_that("set_screen_type honors requested multiplex type when checks favor fixed-pair", {
  screen <- GIScores(
    make_fixed_pair_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = FALSE,
    query_sufficient = FALSE,
    library_sufficient = FALSE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = FALSE
  )
  screen@metadata$requested_screen_type <- "multiplex"

  expect_warning(
    result <- set_screen_type(as(screen, "ScreenBase")),
    "overrides inferred screen type 'fixed_pair'"
  )

  expect_s4_class(result, "MultiplexScreen")
  expect_equal(result@guideGIs@space, c("query_gene", "library_gene"))
  expect_identical(result@metadata$inferred_screen_type, "fixed_pair")
  expect_identical(result@metadata$selected_screen_type, "multiplex")
})

test_that("set_screen_type restores full replicate metadata from fixed-pair style input", {
  screen <- GIScores(
    make_fixed_pair_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = FALSE,
    stable_library_size = FALSE,
    sufficient_tests_per_query = FALSE
  )

  result <- set_screen_type(as(screen, "ScreenBase"))

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
  expect_equal(
    result@guideGIs@replicates,
    c("guide_pair", "tech_rep", "bio_rep")
  )
})
