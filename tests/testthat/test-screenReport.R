make_multiplex_like_scores_for_screenReport <- function() {
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

make_screen_for_screenReport <- function(
  class = c("MultiplexScreen", "FixedPairScreen", "PosAgnMultiplexScreen"),
  checks = list(
    gene_sets_equal = TRUE,
    query_sufficient = TRUE,
    library_sufficient = FALSE,
    stable_library_size = TRUE,
    sufficient_tests_per_query = TRUE,
    avg_tests_per_query = 23
  ),
  failed_queries = character(),
  gi_errors = list(),
  force_fixed_pair = FALSE,
  gene_gis_length = 0L,
  replicate_layers = c("guide_pair", "tech_rep"),
  block_layer = "guide_pair"
) {
  construction_block_layer <- if (length(block_layer) == 0) {
    "guide_pair"
  } else {
    block_layer
  }

  base_screen <- GIScores(
    make_multiplex_like_scores_for_screenReport(),
    block_layer = construction_block_layer
  )

  screen <- as(base_screen, match.arg(class))
  screen@checks <- checks
  screen@errors$query_genes_not_usable <- failed_queries
  screen@errors$GI_computation_errors <- gi_errors
  screen@metadata$force_fixed_pair <- force_fixed_pair
  screen@guideGIs@replicates <- replicate_layers
  screen@guideGIs@block_layer <- block_layer
  screen@geneGIs <- array(numeric(gene_gis_length), dim = gene_gis_length)

  if (identical(class(screen)[1], "PosAgnMultiplexScreen")) {
    screen@symmGeneGIs <- data.table::data.table(
      gene_pair = character(),
      query_gene = character(),
      library_gene = character(),
      GI = numeric(),
      GI_z = numeric(),
      pval = numeric(),
      FDR = numeric()
    )
  }

  screen
}

test_that("screenReport returns structured report content in interactive mode", {
  screen <- make_screen_for_screenReport(
    class = "MultiplexScreen",
    gene_gis_length = 3L,
    gi_errors = list(G1 = NULL, G2 = simpleError("boom"))
  )

  result <- screenReport(screen, interactive = TRUE, print = FALSE)

  expect_type(result, "list")
  expect_named(result, c("overview", "decisions", "checks", "problems"))
  expect_true(any(grepl(
    "Screen class: MultiplexScreen",
    result$overview,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Interpreted design: multiplex",
    result$overview,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Forced fixed-pair run: FALSE",
    result$decisions,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Model objects available:",
    result$decisions,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Gene sets sufficiently overlapping: OK",
    result$checks,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Library size sufficient: PROBLEM",
    result$checks,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Stored model errors: 1",
    result$problems,
    fixed = TRUE
  )))
})

test_that("screenReport prints section headings and key summary lines", {
  screen <- make_screen_for_screenReport(class = "FixedPairScreen")

  output <- paste(
    capture.output(screenReport(screen, print = TRUE)),
    collapse = "\n"
  )

  expect_match(output, "CeRberus screen report")
  expect_match(output, "OVERVIEW")
  expect_match(output, "DECISIONS")
  expect_match(output, "CHECKS")
  expect_match(output, "PROBLEMS")
  expect_match(output, "Screen class: FixedPairScreen")
  expect_match(output, "Selected model strategy: fixed-pair")
})

test_that("screenReport uses fallback text when checks and metadata fields are missing", {
  screen <- make_screen_for_screenReport(
    class = "MultiplexScreen",
    checks = list(),
    failed_queries = NULL,
    gi_errors = list(),
    gene_gis_length = 0L,
    replicate_layers = character(),
    block_layer = character()
  )

  screen@screen_attr$n_query_genes <- NA_real_
  screen@screen_attr$n_lib_genes <- numeric()
  screen@screen_attr$n_all_genes <- NULL
  screen@screen_attr$query_genes_not_in_lib <- character()
  screen@screen_attr$library_genes_not_in_query <- character()
  screen@screen_attr$all_pairs <- character()
  screen@screen_attr$unique_pairs <- character()

  result <- screenReport(screen, interactive = TRUE, print = FALSE)

  expect_equal(result$checks, "No screen checks stored.")
  expect_true(any(grepl(
    "Query genes: not available",
    result$overview,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Library genes: not available",
    result$overview,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "All genes: not available",
    result$overview,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Guide-level replicate layers: not available",
    result$overview,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Duplicate-correlation block layer: not available",
    result$overview,
    fixed = TRUE
  )))
})

test_that("screenReport summarizes failed queries and truncates long affected-query lists", {
  failed_queries <- paste0("Q", seq_len(10L))
  screen <- make_screen_for_screenReport(
    class = "PosAgnMultiplexScreen",
    failed_queries = failed_queries,
    gi_errors = purrr::set_names(
      c(
        list(simpleError("bad1")),
        rep(list(NULL), 8L),
        list(simpleError("bad2"))
      ),
      failed_queries
    ),
    force_fixed_pair = TRUE
  )

  result <- screenReport(screen, interactive = TRUE, print = FALSE)

  expect_true(any(grepl(
    "Position-agnostic output: TRUE",
    result$decisions,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Query genes without usable model: 10",
    result$problems,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Stored model errors: 2",
    result$problems,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Affected query genes: Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8",
    result$problems,
    fixed = TRUE
  )))
  expect_true(any(grepl("+2 more", result$problems, fixed = TRUE)))
})
