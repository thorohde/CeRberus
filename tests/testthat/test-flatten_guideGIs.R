make_gRNA_GI_for_flatten <- function(
  data,
  space = "gene_pair",
  replicates = c("guide_pair", "bio_rep", "tech_rep"),
  block_layer = character(),
  blocks = "none"
) {
  methods::new(
    "gRNA_GI",
    data = data,
    space = space,
    replicates = replicates,
    block_layer = block_layer,
    blocks = blocks,
    use_blocks = FALSE,
    block_description = character(),
    collapse = character()
  )
}

make_fixed_pair_flatten_array <- function() {
  array(
    seq_len(16L),
    dim = c(2L, 2L, 2L, 2L),
    dimnames = list(
      c("A;C", "B;D"),
      c("g1", "g2"),
      c("b1", "b2"),
      c("t1", "t2")
    )
  )
}

make_multiplex_flatten_array <- function() {
  array(
    seq_len(24L),
    dim = c(2L, 3L, 2L, 2L),
    dimnames = list(
      c("Q1", "Q2"),
      c("L1", "L2", "L3"),
      c("g1", "g2"),
      c("b1", "b2")
    )
  )
}

test_that("flatten_guide_gis flattens fixed-pair replicate dimensions into columns", {
  data <- make_fixed_pair_flatten_array()
  object <- make_gRNA_GI_for_flatten(data)

  result <- flatten_guide_gis(object)

  expect_s4_class(result, "gRNA_GI")
  expect_equal(dim(result@data), c(2L, 8L))
  expect_equal(rownames(result@data), c("A;C", "B;D"))
  expect_equal(
    colnames(result@data),
    c(
      "g1_b1_t1",
      "g1_b1_t2",
      "g1_b2_t1",
      "g1_b2_t2",
      "g2_b1_t1",
      "g2_b1_t2",
      "g2_b2_t1",
      "g2_b2_t2"
    )
  )
  expect_equal(result@data["A;C", "g1_b1_t1"], data["A;C", "g1", "b1", "t1"])
  expect_equal(result@data["B;D", "g2_b2_t2"], data["B;D", "g2", "b2", "t2"])
})

test_that("flatten_guide_gis flattens multiplex replicate dimensions while preserving two space axes", {
  data <- make_multiplex_flatten_array()
  object <- make_gRNA_GI_for_flatten(
    data,
    space = c("query_gene", "library_gene"),
    replicates = c("guide_pair", "bio_rep")
  )

  result <- flatten_guide_gis(object)

  expect_equal(dim(result@data), c(2L, 3L, 4L))
  expect_equal(dimnames(result@data)[[1L]], c("Q1", "Q2"))
  expect_equal(dimnames(result@data)[[2L]], c("L1", "L2", "L3"))
  expect_equal(dimnames(result@data)[[3L]], c("g1_b1", "g1_b2", "g2_b1", "g2_b2"))
  expect_equal(result@data["Q1", "L1", "g1_b1"], data["Q1", "L1", "g1", "b1"])
  expect_equal(result@data["Q2", "L3", "g2_b2"], data["Q2", "L3", "g2", "b2"])
})

test_that("flatten_guide_gis derives or preserves block labels", {
  data <- make_fixed_pair_flatten_array()
  descriptions <- c(
    "g1_b1_t1", "g1_b1_t2", "g1_b2_t1", "g1_b2_t2",
    "g2_b1_t1", "g2_b1_t2", "g2_b2_t1", "g2_b2_t2"
  )
  cases <- list(
    guide_pair = list(
      object = make_gRNA_GI_for_flatten(
        data,
        blocks = character(),
        block_layer = "guide_pair"
      ),
      blocks = c("g1", "g1", "g1", "g1", "g2", "g2", "g2", "g2")
    ),
    bio_rep = list(
      object = make_gRNA_GI_for_flatten(
        data,
        blocks = character(),
        block_layer = "bio_rep"
      ),
      blocks = c("b1", "b1", "b2", "b2", "b1", "b1", "b2", "b2")
    ),
    tech_rep = list(
      object = make_gRNA_GI_for_flatten(
        data,
        blocks = character(),
        block_layer = "tech_rep"
      ),
      blocks = c("t1", "t2", "t1", "t2", "t1", "t2", "t1", "t2")
    ),
    explicit_none = list(
      object = make_gRNA_GI_for_flatten(
        data,
        block_layer = "guide_pair",
        blocks = "none"
      ),
      blocks = "none"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    result <- flatten_guide_gis(case$object)
    expect_equal(
      list(
        use_blocks = result@use_blocks,
        description = result@block_description,
        blocks = result@blocks
      ),
      list(use_blocks = TRUE, description = descriptions, blocks = case$blocks),
      info = case_name
    )
  })
})

test_that("flatten_guide_gis does not use blocks when block layer is missing or sole replicate", {
  no_block_layer <- make_gRNA_GI_for_flatten(
    make_fixed_pair_flatten_array(),
    block_layer = character()
  )
  single_replicate <- make_gRNA_GI_for_flatten(
    array(
      seq_len(4L),
      dim = c(2L, 2L),
      dimnames = list(c("A;C", "B;D"), c("g1", "g2"))
    ),
    replicates = "guide_pair",
    block_layer = "guide_pair"
  )

  no_block_result <- flatten_guide_gis(no_block_layer)
  single_replicate_result <- flatten_guide_gis(single_replicate)

  expect_false(isTRUE(no_block_result@use_blocks))
  expect_false(single_replicate_result@use_blocks)
})
