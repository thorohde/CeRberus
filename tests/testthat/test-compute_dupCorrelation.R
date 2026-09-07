make_gRNA_GI_for_dupCorrelation <- function(
  data,
  space = "gene_pair",
  replicates = "replicate",
  blocks = character(),
  use_blocks = FALSE
) {
  methods::new(
    "gRNA_GI",
    data = data,
    space = space,
    replicates = replicates,
    block_layer = character(),
    blocks = blocks,
    use_blocks = use_blocks,
    block_description = character(),
    collapse = character()
  )
}

make_screen_for_dupCorrelation <- function(guideGIs) {
  methods::new(
    "ScreenBase",
    guideLFCs = methods::new(
      "gRNA_LFC",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character()
    ),
    guideGIs = guideGIs,
    limma_models = list(),
    geneGIs = array(numeric(), dim = 0),
    screen_attr = methods::new("ScreenDesign"),
    dupCorrelation = numeric(),
    metadata = list(),
    checks = list(),
    errors = list()
  )
}

with_mocked_duplicateCorrelation <- function(code, correlations) {
  calls <- new.env(parent = emptyenv())
  calls$args <- list()
  calls$correlations <- correlations

  testthat::local_mocked_bindings(
    duplicateCorrelation = function(object, block, ndups, ...) {
      calls$args[[length(calls$args) + 1L]] <- list(
        object = object,
        block = if (missing(block)) NULL else block,
        ndups = if (missing(ndups)) NULL else ndups
      )

      list(
        consensus.correlation = calls$correlations[[length(calls$args)]]
      )
    },
    .package = "limma"
  )

  result <- force(code)
  list(result = result, calls = calls$args)
}

make_fixed_pair_dupCorrelation_matrix <- function() {
  matrix(
    seq_len(12L),
    nrow = 3L,
    dimnames = list(
      gene_pair = c("A;C", "B;D", "E;F"),
      replicate = paste0("r", seq_len(4L))
    )
  )
}

make_multiplex_dupCorrelation_array <- function() {
  array(
    seq_len(24L),
    dim = c(2L, 3L, 4L),
    dimnames = list(
      query_gene = c("Q1", "Q2"),
      library_gene = c("L1", "L2", "L3"),
      replicate = paste0("r", seq_len(4L))
    )
  )
}


make_global_dupCorrelation_scores <- function(n_genes = 20L) {
  genes <- paste0("G", seq_len(n_genes))
  input <- expand.grid(
    query_gene = genes,
    library_gene = genes,
    guide_pair = c("g1", "g2"),
    bio_rep = c("b1", "b2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  input$GI <- seq_len(nrow(input))
  input
}

test_that("compute_dup_correlation handles fixed-pair blocking variants", {
  data <- make_fixed_pair_dupCorrelation_matrix()
  blocks <- c("b1", "b1", "b2", "b2")
  cases <- list(
    unblocked = list(
      use_blocks = FALSE,
      correlation = 0.123,
      block = NULL,
      ndups = 1
    ),
    blocked = list(
      use_blocks = TRUE,
      correlation = 0.456,
      block = blocks,
      ndups = NULL
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    object <- make_gRNA_GI_for_dupCorrelation(
      data,
      blocks = if (case$use_blocks) blocks else character(),
      use_blocks = case$use_blocks
    )
    mocked <- with_mocked_duplicateCorrelation(
      compute_dup_correlation(object),
      correlations = list(case$correlation)
    )

    expect_equal(mocked$result, case$correlation, info = case_name)
    expect_equal(mocked$calls, list(list(
      object = data,
      block = case$block,
      ndups = case$ndups
    )), info = case_name)
  })
})

test_that("compute_dup_correlation handles multiplex blocking variants", {
  data <- make_multiplex_dupCorrelation_array()
  blocks <- c("b1", "b1", "b2", "b2")
  cases <- list(
    unblocked = list(
      use_blocks = FALSE,
      correlations = list(0.1, 0.2),
      result = c(Q1 = 0.1, Q2 = 0.2),
      block = NULL,
      ndups = 1
    ),
    blocked = list(
      use_blocks = TRUE,
      correlations = list(0.3, 0.4),
      result = c(Q1 = 0.3, Q2 = 0.4),
      block = blocks,
      ndups = NULL
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    object <- make_gRNA_GI_for_dupCorrelation(
      data,
      space = c("query_gene", "library_gene"),
      replicates = "replicate",
      blocks = if (case$use_blocks) blocks else character(),
      use_blocks = case$use_blocks
    )
    mocked <- with_mocked_duplicateCorrelation(
      compute_dup_correlation(object),
      correlations = case$correlations
    )
    expected_calls <- purrr::map(c("Q1", "Q2"), function(query) {
      list(
        object = data[query, , ],
        block = case$block,
        ndups = case$ndups
      )
    })

    expect_equal(mocked$result, case$result, info = case_name)
    expect_equal(mocked$calls, expected_calls, info = case_name)
  })
})

test_that("compute_dup_correlation for ScreenBase stores the guide-level duplicate correlation", {
  data <- make_multiplex_dupCorrelation_array()
  guideGIs <- make_gRNA_GI_for_dupCorrelation(
    data,
    space = c("query_gene", "library_gene"),
    replicates = "replicate"
  )
  screen <- make_screen_for_dupCorrelation(guideGIs)

  mocked <- with_mocked_duplicateCorrelation(
    compute_dup_correlation(screen),
    correlations = list(0.5, 0.6)
  )

  expect_s4_class(mocked$result, "ScreenBase")
  expect_equal(mocked$result@dupCorrelation, c(Q1 = 0.5, Q2 = 0.6))
  expect_equal(mocked$result@guideGIs, guideGIs)
  expect_length(mocked$calls, 2L)
})


test_that("global_preaverage computes one duplicate correlation over all pairs", {
  screen <- GIScores(
    make_global_dupCorrelation_scores(),
    pos_agnostic = TRUE,
    symmetric_analysis_method = "global_preaverage",
    block_layer = "guide_pair"
  )
  expected_data <- screen@guideGIs@data
  expected_blocks <- screen@guideGIs@blocks

  mocked <- with_mocked_duplicateCorrelation(
    compute_dup_correlation(screen),
    correlations = list(0.321)
  )

  expect_s4_class(mocked$result, "PosAgnMultiplexScreen")
  expect_equal(mocked$result@dupCorrelation, 0.321)
  expect_length(mocked$calls, 1L)
  expect_equal(mocked$calls[[1L]]$object, expected_data)
  expect_equal(mocked$calls[[1L]]$block, expected_blocks)
  expect_null(mocked$calls[[1L]]$ndups)
})
