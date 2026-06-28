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

test_that("gRNA_GI stores provided data and metadata slots", {
  data <- array(
    seq_len(8L),
    dim = c(2L, 2L, 2L),
    dimnames = list(
      query_gene = c("A", "B"),
      library_gene = c("C", "D"),
      guide_pair = c("g1", "g2")
    )
  )

  result <- make_gRNA_GI_for_class_tests(
    data = data,
    space = c("query_gene", "library_gene"),
    replicates = "guide_pair",
    block_layer = "guide_pair",
    blocks = c("g1", "g2"),
    use_blocks = TRUE,
    block_description = c("g1", "g2"),
    collapse = character()
  )

  expect_s4_class(result, "gRNA_GI")
  expect_equal(result@data, data)
  expect_equal(result@space, c("query_gene", "library_gene"))
  expect_equal(result@replicates, "guide_pair")
  expect_equal(result@block_layer, "guide_pair")
  expect_equal(result@blocks, c("g1", "g2"))
  expect_true(result@use_blocks)
  expect_equal(result@block_description, c("g1", "g2"))
  expect_equal(result@collapse, character())
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
  expect_length(result@data, 0L)
  expect_equal(result@space, character())
  expect_equal(result@replicates, character())
  expect_equal(result@block_layer, character())
  expect_equal(result@blocks, character())
  expect_false(result@use_blocks)
  expect_equal(result@block_description, character())
  expect_equal(result@collapse, character())
})

test_that("gRNA_GI preserves array dimnames and dimension order", {
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

  expect_equal(unname(dim(result@data)), c(2L, 3L, 2L))
  expect_equal(
    names(dimnames(result@data)),
    c("query_gene", "library_gene", "tech_rep")
  )
  expect_equal(dimnames(result@data)$query_gene, c("Q1", "Q2"))
  expect_equal(dimnames(result@data)$library_gene, c("L1", "L2", "L3"))
  expect_equal(dimnames(result@data)$tech_rep, c("t1", "t2"))
})

test_that("gRNA_GI validates declared slot types", {
  expect_error(
    methods::new(
      "gRNA_GI",
      data = 1,
      space = character(),
      replicates = character(),
      block_layer = character(),
      blocks = character(),
      use_blocks = FALSE,
      block_description = character(),
      collapse = character()
    )
  )

  expect_error(
    methods::new(
      "gRNA_GI",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character(),
      block_layer = character(),
      blocks = character(),
      use_blocks = "FALSE",
      block_description = character(),
      collapse = character()
    )
  )
})
