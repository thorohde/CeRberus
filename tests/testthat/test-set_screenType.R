make_multiplex_like_scores_for_set_screenType <- function() {
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

make_fixed_pair_like_scores_for_set_screenType <- function() {
  data.frame(
    query_gene = rep(c("A", "B"), each = 4),
    library_gene = rep(c("C", "D"), each = 4),
    bio_rep = rep(c("b1", "b2"), 4),
    tech_rep = rep(c("t1", "t2"), each = 2, times = 2),
    guide_pair = rep(c("g1", "g2"), 4),
    GI = seq_len(8)
  )
}

test_that("set_screenType keeps multiplex-compatible screens as multiplex", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screenType(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = TRUE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE
  )

  result <- set_screenType(as(screen, "ScreenBase"))

  expect_s4_class(result, "MultiplexScreen")
  expect_equal(result@guideGIs@space, c("query_gene", "library_gene"))
  expect_equal(result@guideGIs@replicates, "guide_pair")
  expect_equal(result@metadata$input, screen@metadata$input)
})

test_that("set_screenType switches to fixed-pair when library checks fail", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screenType(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = FALSE,
    stable_library_size = FALSE,
    sufficient_tests_per_query = FALSE
  )

  result <- set_screenType(as(screen, "ScreenBase"))

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
  expect_equal(result@guideGIs@replicates, "guide_pair")
})

test_that("set_screenType warns and falls back to fixed-pair for unknown designs", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screenType(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = FALSE,
    library_sufficient = TRUE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE
  )

  expect_warning(
    result <- set_screenType(as(screen, "ScreenBase")),
    "Unknown screen design! Forcing fixed pair run"
  )

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
})

test_that("set_screenType honors force_fixed_pair even when checks allow multiplex", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screenType(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = TRUE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE
  )
  screen@metadata$force_fixed_pair <- TRUE

  expect_warning(
    result <- set_screenType(as(screen, "ScreenBase")),
    "Set up to use fixed pair structure"
  )

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
  expect_equal(result@guideGIs@replicates, "guide_pair")
})

test_that("set_screenType restores full replicate metadata from fixed-pair style input", {
  screen <- GIScores(
    make_fixed_pair_like_scores_for_set_screenType(),
    block_layer = "guide_pair"
  )

  screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = FALSE,
    stable_library_size = FALSE,
    sufficient_tests_per_query = FALSE
  )

  result <- set_screenType(as(screen, "ScreenBase"))

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
  expect_equal(
    result@guideGIs@replicates,
    c("guide_pair", "tech_rep", "bio_rep")
  )
})
