make_screen_for_create_log <- function() {
  methods::new(
    "MultiplexScreen",
    guideLFCs = methods::new(
      "gRNA_LFC",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character()
    ),
    guideGIs = methods::new(
      "gRNA_GI",
      data = array(
        c(1, NA_real_, 3, 4),
        dim = c(2, 2),
        dimnames = list(
          query_gene = c("Q1", "Q2"),
          replicate = c("g1", "g2")
        )
      ),
      space = "query_gene",
      replicates = c("guide_pair", "bio_rep"),
      block_layer = "guide_pair",
      blocks = c("g1", "g2"),
      use_blocks = TRUE,
      block_description = c("g1", "g2"),
      collapse = "tech_rep"
    ),
    limma_models = list(Q1 = list(coefficients = 1), Q2 = NULL),
    geneGIs = array(numeric(), dim = 0),
    screen_attr = list(
      n_query_genes = 2L,
      n_lib_genes = 3L,
      n_all_genes = 5L,
      all_pairs = c("Q1;L1", "Q1;L2", "Q2;L3"),
      unique_pairs = c("L1;Q1", "L2;Q1", "L3;Q2")
    ),
    dupCorrelation = c(Q1 = 0.1, Q2 = NA_real_),
    metadata = list(
      config = "guide_pair_used",
      input = data.table::data.table(
        query_gene = "SECRET_GENE",
        library_gene = "SECRET_LIBRARY",
        GI = 999
      )
    ),
    checks = list(
      gene_sets_equal = FALSE,
      sufficient_tests_per_query = TRUE
    ),
    errors = list(
      query_genes_not_usable = c("Q2", "Q3", "Q4"),
      GI_computation_errors = list(
        Q1 = NULL,
        Q2 = simpleError("model failed for Q2")
      )
    )
  )
}

test_that("create_log records bounded diagnostics without raw input values", {
  result <- create_log(
    make_screen_for_create_log(),
    status = "failed",
    stage = "fit_models",
    condition = simpleError("pipeline exploded"),
    max_items = 2L
  )

  expect_type(result, "character")
  expect_length(result, 1L)
  expect_match(result, "Status: failed", fixed = TRUE)
  expect_match(result, "Pipeline stage: fit_models", fixed = TRUE)
  expect_match(result, "Condition: pipeline exploded", fixed = TRUE)
  expect_match(result, "Configuration: guide_pair_used", fixed = TRUE)
  expect_match(result, "Class: MultiplexScreen", fixed = TRUE)
  expect_match(
    result,
    "Guide GI data: dimensions=2 x 2; values=4; missing=1",
    fixed = TRUE
  )
  expect_match(
    result,
    "Duplicate correlation: values=2; finite=1; missing=1",
    fixed = TRUE
  )
  expect_match(
    result,
    "Failed queries (3): Q2 | Q3 ... +1 more",
    fixed = TRUE
  )
  expect_match(result, "model failed for Q2", fixed = TRUE)
  expect_false(grepl("SECRET_GENE|SECRET_LIBRARY|999", result))
})

test_that("create_log handles an empty partially constructed screen", {
  result <- create_log(methods::new("ScreenBase"))

  expect_match(result, "Status: available", fixed = TRUE)
  expect_match(result, "Guide GI data: empty", fixed = TRUE)
  expect_match(result, "Duplicate correlation: not available", fixed = TRUE)
  expect_match(result, "Stored model errors (0): none", fixed = TRUE)
})

test_that("create_log validates controls", {
  screen <- methods::new("ScreenBase")

  expect_error(create_log(screen, status = NA_character_), "status must")
  expect_error(create_log(screen, stage = character()), "stage must")
  expect_error(create_log(screen, max_items = 0L), "max_items must")
})
