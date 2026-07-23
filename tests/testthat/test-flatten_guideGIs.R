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

test_that("flatten_guide_gis extracts block labels from guide-pair replicate descriptions", {
  object <- make_gRNA_GI_for_flatten(
    make_fixed_pair_flatten_array(),
    blocks = character(),
    block_layer = "guide_pair"
  )

  result <- flatten_guide_gis(object)

  expect_true(result@use_blocks)
  expect_equal(result@block_description, colnames(result@data))
  expect_equal(result@blocks, c("g1", "g1", "g1", "g1", "g2", "g2", "g2", "g2"))
})

test_that("flatten_guide_gis extracts block labels from bio-rep and tech-rep descriptions", {
  bio_object <- make_gRNA_GI_for_flatten(
    make_fixed_pair_flatten_array(),
    blocks = character(),
    block_layer = "bio_rep"
  )
  tech_object <- make_gRNA_GI_for_flatten(
    make_fixed_pair_flatten_array(),
    blocks = character(),
    block_layer = "tech_rep"
  )

  bio_result <- flatten_guide_gis(bio_object)
  tech_result <- flatten_guide_gis(tech_object)

  expect_equal(bio_result@blocks, c("b1", "b1", "b2", "b2", "b1", "b1", "b2", "b2"))
  expect_equal(tech_result@blocks, c("t1", "t2", "t1", "t2", "t1", "t2", "t1", "t2"))
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

test_that("flatten_guide_gis preserves explicit single-value blocks such as none", {
  object <- make_gRNA_GI_for_flatten(
    make_fixed_pair_flatten_array(),
    block_layer = "guide_pair",
    blocks = "none"
  )

  result <- flatten_guide_gis(object)

  expect_true(result@use_blocks)
  expect_equal(result@blocks, "none")
})
