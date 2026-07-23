make_multiplex_like_scores_for_screen_report <- function() {
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

make_screen_for_screen_report <- function(
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
  requested_screen_type = "auto",
  inferred_screen_type = NULL,
  selected_screen_type = NULL,
  gene_gis_length = 0L,
  replicate_layers = c("guide_pair", "tech_rep"),
  block_layer = "guide_pair",
  metadata = list(),
  limma_models = list()
) {
  construction_block_layer <- if (length(block_layer) == 0) {
    "guide_pair"
  } else {
    block_layer
  }

  base_screen <- GIScores(
    make_multiplex_like_scores_for_screen_report(),
    block_layer = construction_block_layer
  )

  screen <- as(base_screen, match.arg(class))
  class_screen_type <- switch(
    class(screen)[1L],
    FixedPairScreen = "fixed_pair",
    MultiplexScreen = "multiplex",
    PosAgnMultiplexScreen = "multiplex"
  )
  if (is.null(inferred_screen_type)) {
    inferred_screen_type <- class_screen_type
  }
  if (is.null(selected_screen_type)) {
    selected_screen_type <- class_screen_type
  }
  screen@checks <- checks
  screen@errors$query_genes_not_usable <- failed_queries
  screen@errors$GI_computation_errors <- gi_errors
  screen@metadata$requested_screen_type <- requested_screen_type
  screen@metadata$inferred_screen_type <- inferred_screen_type
  screen@metadata$selected_screen_type <- selected_screen_type
  screen@metadata <- utils::modifyList(screen@metadata, metadata)
  screen@limma_models <- limma_models
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

test_that("screen_report returns structured report content in interactive mode", {
  screen <- make_screen_for_screen_report(
    class = "MultiplexScreen",
    gene_gis_length = 3L,
    gi_errors = list(G1 = NULL, G2 = simpleError("boom"))
  )

  result <- screen_report(screen, interactive = TRUE, print = FALSE)

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
    "Requested screen type: auto",
    result$decisions,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Inferred screen type: multiplex",
    result$decisions,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Selected screen type: multiplex",
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

test_that("screen_report prints section headings and key summary lines", {
  screen <- make_screen_for_screen_report(class = "FixedPairScreen")

  output <- paste(
    capture.output(screen_report(screen, print = TRUE)),
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

test_that("screen_report uses fallback text when checks and metadata fields are missing", {
  screen <- make_screen_for_screen_report(
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

  result <- screen_report(screen, interactive = TRUE, print = FALSE)

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

test_that("screen_report summarizes failed queries and truncates long affected-query lists", {
  failed_queries <- paste0("Q", seq_len(10L))
  screen <- make_screen_for_screen_report(
    class = "PosAgnMultiplexScreen",
    failed_queries = failed_queries,
    gi_errors = purrr::set_names(
      c(
        list(simpleError("bad1")),
        rep(list(NULL), 8L),
        list(simpleError("bad2"))
      ),
      failed_queries
    )
  )

  result <- screen_report(screen, interactive = TRUE, print = FALSE)

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


test_that("screen_report identifies one global position-agnostic model", {
  fit <- limma::lmFit(matrix(c(1, 2, 3, 4), nrow = 2L))
  screen <- make_screen_for_screen_report(
    class = "PosAgnMultiplexScreen",
    metadata = list(symmetric_analysis_method = "global_preaverage"),
    limma_models = fit
  )

  result <- screen_report(screen, interactive = TRUE, print = FALSE)

  expect_true(any(grepl(
    "Symmetric analysis method: global_preaverage",
    result$decisions,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Model objects available: 1",
    result$decisions,
    fixed = TRUE
  )))
})
