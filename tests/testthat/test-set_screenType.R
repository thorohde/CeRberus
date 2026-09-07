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

test_that("set_screen_type infers screen type from library checks", {
  screen <- GIScores(
    make_multiplex_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )
  cases <- list(
    multiplex_compatible = list(
      checks = list(
        gene_sets_equal = TRUE,
        query_sufficient = TRUE,
        library_sufficient = TRUE,
        stable_library_size = TRUE,
        sufficient_tests_per_query = TRUE
      ),
      class = "MultiplexScreen",
      space = c("query_gene", "library_gene")
    ),
    library_checks_fail = list(
      checks = list(
        gene_sets_equal = TRUE,
        query_sufficient = TRUE,
        library_sufficient = FALSE,
        stable_library_size = FALSE,
        sufficient_tests_per_query = FALSE
      ),
      class = "FixedPairScreen",
      space = "gene_pair"
    ),
    query_checks_fail = list(
      checks = list(
        gene_sets_equal = TRUE,
        query_sufficient = FALSE,
        library_sufficient = TRUE,
        stable_library_size = TRUE,
        sufficient_tests_per_query = TRUE
      ),
      class = "FixedPairScreen",
      space = "gene_pair"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    candidate <- screen
    candidate@checks <- case$checks
    result <- set_screen_type(as(candidate, "ScreenBase"))

    expect_equal(
      list(
        class = class(result)[[1]],
        space = result@guideGIs@space,
        replicates = result@guideGIs@replicates,
        input = result@metadata$input
      ),
      list(
        class = case$class,
        space = case$space,
        replicates = "guide_pair",
        input = screen@metadata$input
      ),
      info = case_name
    )
  })
})

test_that("set_screen_type honors explicit screen-type overrides", {
  multiplex_screen <- GIScores(
    make_multiplex_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )
  multiplex_screen@checks <- list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = TRUE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE
  )
  fixed_pair_screen <- GIScores(
    make_fixed_pair_like_scores_for_set_screen_type(),
    block_layer = "guide_pair"
  )
  fixed_pair_screen@checks <- list(
    gene_sets_equal = FALSE,
    query_sufficient = FALSE,
    library_sufficient = FALSE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = FALSE
  )
  cases <- list(
    fixed_pair = list(
      screen = multiplex_screen,
      requested = "fixed_pair",
      inferred = "multiplex",
      class = "FixedPairScreen",
      space = "gene_pair"
    ),
    multiplex = list(
      screen = fixed_pair_screen,
      requested = "multiplex",
      inferred = "fixed_pair",
      class = "MultiplexScreen",
      space = c("query_gene", "library_gene")
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    screen <- case$screen
    screen@metadata$requested_screen_type <- case$requested

    expect_warning(
      result <- set_screen_type(as(screen, "ScreenBase")),
      paste0("overrides inferred screen type '", case$inferred, "'"),
      info = case_name
    )
    expect_equal(
      list(
        class = class(result)[[1]],
        space = result@guideGIs@space,
        inferred = result@metadata$inferred_screen_type,
        selected = result@metadata$selected_screen_type
      ),
      list(
        class = case$class,
        space = case$space,
        inferred = case$inferred,
        selected = case$requested
      ),
      info = case_name
    )
  })
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

  expect_equal(
    list(
      class = class(result)[[1]],
      space = result@guideGIs@space,
      replicates = result@guideGIs@replicates
    ),
    list(
      class = "FixedPairScreen",
      space = "gene_pair",
      replicates = c("guide_pair", "tech_rep", "bio_rep")
    )
  )
})
