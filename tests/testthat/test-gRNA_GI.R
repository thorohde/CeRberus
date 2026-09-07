make_gRNA_GI_for_class_tests <- function(
  data = array(
    seq_len(8L),
    dim = c(2L, 2L, 2L),
    dimnames = list(
      gene_pair = c("A;C", "B;D"),
      guide_pair = c("g1", "g2"),
      bio_rep = c("b1", "b2")
    )
  ),
  space = "gene_pair",
  replicates = c("guide_pair", "bio_rep"),
  block_layer = "guide_pair",
  blocks = c("g1", "g2"),
  use_blocks = TRUE,
  block_description = c("g1_b1", "g2_b2"),
  collapse = "bio_rep"
) {
  methods::new(
    "gRNA_GI",
    data = data,
    space = space,
    replicates = replicates,
    block_layer = block_layer,
    blocks = blocks,
    use_blocks = use_blocks,
    block_description = block_description,
    collapse = collapse
  )
}

test_that("gRNA_GI stores valid array and metadata slots", {
  data <- array(
    seq_len(12L),
    dim = c(2L, 3L, 2L),
    dimnames = list(
      query_gene = c("Q1", "Q2"),
      library_gene = c("L1", "L2", "L3"),
      tech_rep = c("t1", "t2")
    )
  )
  result <- make_gRNA_GI_for_class_tests(
    data = data,
    space = c("query_gene", "library_gene"),
    replicates = "tech_rep",
    block_layer = character(),
    blocks = character(),
    use_blocks = FALSE,
    block_description = character(),
    collapse = character()
  )

  expect_s4_class(result, "gRNA_GI")
  expect_equal(
    list(
      data = result@data,
      space = result@space,
      replicates = result@replicates,
      block_layer = result@block_layer,
      blocks = result@blocks,
      use_blocks = result@use_blocks,
      block_description = result@block_description,
      collapse = result@collapse
    ),
    list(
      data = data,
      space = c("query_gene", "library_gene"),
      replicates = "tech_rep",
      block_layer = character(),
      blocks = character(),
      use_blocks = FALSE,
      block_description = character(),
      collapse = character()
    )
  )
})

test_that("gRNA_GI supports empty defaults for optional metadata slots", {
  result <- methods::new(
    "gRNA_GI",
    data = array(numeric(), dim = 0),
    space = character(),
    replicates = character(),
    block_layer = character(),
    blocks = character(),
    use_blocks = FALSE,
    block_description = character(),
    collapse = character()
  )

  expect_s4_class(result, "gRNA_GI")
  expect_equal(
    list(
      data = result@data,
      space = result@space,
      replicates = result@replicates,
      block_layer = result@block_layer,
      blocks = result@blocks,
      use_blocks = result@use_blocks,
      block_description = result@block_description,
      collapse = result@collapse
    ),
    list(
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character(),
      block_layer = character(),
      blocks = character(),
      use_blocks = FALSE,
      block_description = character(),
      collapse = character()
    )
  )
})

test_that("gRNA_GI rejects invalid slot and array metadata", {
  cases <- list(
    non_array_data = list(
      call = function() methods::new(
        "gRNA_GI",
        data = 1,
        space = character(),
        replicates = character(),
        block_layer = character(),
        blocks = character(),
        use_blocks = FALSE,
        block_description = character(),
        collapse = character()
      ),
      error = NULL
    ),
    non_logical_blocks = list(
      call = function() methods::new(
        "gRNA_GI",
        data = array(numeric(), dim = 0),
        space = character(),
        replicates = character(),
        block_layer = character(),
        blocks = character(),
        use_blocks = "FALSE",
        block_description = character(),
        collapse = character()
      ),
      error = NULL
    ),
    multiple_block_layers = list(
      call = function() make_gRNA_GI_for_class_tests(
        block_layer = c("guide_pair", "bio_rep")
      ),
      error = "block_layer.*empty or one"
    ),
    non_replicate_block = list(
      call = function() make_gRNA_GI_for_class_tests(block_layer = "gene_pair"),
      error = "biological.*block_layer"
    ),
    missing_blocks = list(
      call = function() make_gRNA_GI_for_class_tests(
        blocks = character(),
        use_blocks = TRUE
      ),
      error = "blocks.*populated"
    ),
    mismatched_block_description = list(
      call = function() make_gRNA_GI_for_class_tests(block_description = "one_column"),
      error = "final 'data' dimension"
    ),
    mismatched_blocks = list(
      call = function() make_gRNA_GI_for_class_tests(blocks = "g1"),
      error = "equal lengths"
    ),
    wrong_array_rank = list(
      call = function() make_gRNA_GI_for_class_tests(
        data = array(seq_len(4L), dim = c(2L, 2L)),
        space = "gene_pair",
        replicates = character(),
        block_layer = character(),
        blocks = character(),
        use_blocks = FALSE,
        block_description = character(),
        collapse = character()
      ),
      error = "array rank"
    ),
    duplicate_dimension_names = list(
      call = function() make_gRNA_GI_for_class_tests(
        space = "gene_pair",
        replicates = "gene_pair"
      ),
      error = "unique dimension names"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    if (is.null(case$error)) {
      expect_error(case$call(), info = case_name)
    } else {
      expect_error(case$call(), case$error, info = case_name)
    }
  })
})
