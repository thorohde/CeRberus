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
    screen_attr = make_screen_design(
      query_genes = c("Q1", "Q2"),
      library_genes = c("L1", "L2", "L3"),
      all_pairs = c("Q1;L1", "Q1;L2", "Q2;L3")
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

test_that("create_log records bounded diagnostics without leaking raw input", {
  populated <- create_log(
    make_screen_for_create_log(),
    status = "failed",
    stage = "fit_models",
    condition = simpleError("pipeline exploded"),
    max_items = 2L
  )
  populated_markers <- c(
    "Status: failed",
    "Pipeline stage: fit_models",
    "Condition: pipeline exploded",
    "Configuration: guide_pair_used",
    "Class: MultiplexScreen",
    "Guide GI data: dimensions=2 x 2; values=4; missing=1",
    "Duplicate correlation: values=2; finite=1; missing=1",
    "Failed queries (3): Q2 | Q3 ... +1 more",
    "model failed for Q2"
  )

  expect_true(
    is.character(populated) &&
      length(populated) == 1L &&
      all(vapply(populated_markers, grepl, logical(1), x = populated, fixed = TRUE)) &&
      !grepl("SECRET_GENE|SECRET_LIBRARY|999", populated)
  )

  empty <- create_log(methods::new("ScreenBase"))
  empty_markers <- c(
    "Status: available",
    "Guide GI data: empty",
    "Duplicate correlation: not available",
    "Stored model errors (0): none"
  )

  expect_true(all(vapply(
    empty_markers,
    grepl,
    logical(1),
    x = empty,
    fixed = TRUE
  )))
})

test_that("create_log validates controls", {
  screen <- methods::new("ScreenBase")

  expect_error(
    create_log(screen, status = NA_character_),
    "^status must be a single non-missing character value\\.$"
  )
  expect_error(
    create_log(screen, stage = character()),
    "^stage must be a single non-missing character value\\.$"
  )
  expect_error(
    create_log(screen, max_items = 0L),
    "^max_items must be a single positive whole number\\.$"
  )
})
