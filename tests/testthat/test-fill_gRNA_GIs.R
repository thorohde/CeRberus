make_gRNA_GI_for_fill <- function(
  space = "gene_pair",
  replicates = c("guide_pair", "bio_rep", "tech_rep"),
  block_layer = character(),
  blocks = character(),
  use_blocks = FALSE,
  block_description = character(),
  collapse = character()
) {
  methods::new(
    "gRNA_GI",
    data = array(numeric(), dim = 0),
    space = space,
    replicates = replicates,
    block_layer = block_layer,
    blocks = blocks,
    use_blocks = use_blocks,
    block_description = block_description,
    collapse = collapse
  )
}

make_fixed_pair_fill_scores <- function() {
  expand.grid(
    gene_pair = c("A;C", "B;D"),
    guide_pair = c("g1", "g2"),
    bio_rep = c("b1", "b2"),
    tech_rep = c("t1", "t2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    transform(
      GI = seq_len(16L),
      alternative_score = seq_len(16L) + 100L
    )
}

make_multiplex_fill_scores <- function() {
  expand.grid(
    query_gene = c("G1", "G2"),
    library_gene = c("G1", "G2", "G3"),
    guide_pair = c("g1", "g2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    transform(GI = seq_len(12L))
}

test_that("fill_gRNA_GIs fills fixed-pair guide GI arrays from long input", {
  input <- make_fixed_pair_fill_scores()
  object <- make_gRNA_GI_for_fill()

  result <- CeRberus:::fill_gRNA_GIs(object, input)

  expect_s4_class(result, "gRNA_GI")
  expect_true(is.array(result@data))
  expect_equal(dim(result@data), c(2L, 2L, 2L, 2L))
  expect_equal(dimnames(result@data)[[1L]], c("A;C", "B;D"))
  expect_equal(dimnames(result@data)[[2L]], c("g1", "g2"))
  expect_equal(dimnames(result@data)[[3L]], c("b1", "b2"))
  expect_equal(dimnames(result@data)[[4L]], c("t1", "t2"))

  expect_equal(result@data["A;C", "g1", "b1", "t1"], 1L)
  expect_equal(result@data["B;D", "g2", "b2", "t2"], 16L)
})

test_that("fill_gRNA_GIs fills multiplex guide GI arrays with two space dimensions", {
  input <- make_multiplex_fill_scores()
  object <- make_gRNA_GI_for_fill(
    space = c("query_gene", "library_gene"),
    replicates = "guide_pair"
  )

  result <- CeRberus:::fill_gRNA_GIs(object, input)

  expect_s4_class(result, "gRNA_GI")
  expect_equal(dim(result@data), c(2L, 3L, 2L))
  expect_equal(dimnames(result@data)[[1L]], c("G1", "G2"))
  expect_equal(dimnames(result@data)[[2L]], c("G1", "G2", "G3"))
  expect_equal(dimnames(result@data)[[3L]], c("g1", "g2"))

  expect_equal(result@data["G1", "G1", "g1"], 1L)
  expect_equal(result@data["G2", "G3", "g2"], 12L)
})

test_that("fill_gRNA_GIs supports custom value variables", {
  input <- make_fixed_pair_fill_scores()
  object <- make_gRNA_GI_for_fill(replicates = "guide_pair")
  input <- input[input$bio_rep == "b1" & input$tech_rep == "t1", ]

  result <- CeRberus:::fill_gRNA_GIs(
    object,
    input,
    value_var = "alternative_score"
  )

  expect_equal(result@data["A;C", "g1"], 101L)
  expect_equal(result@data["B;D", "g2"], 104L)
  expect_false(any(result@data == input$GI))
})

test_that("fill_gRNA_GIs only changes the data slot", {
  input <- make_fixed_pair_fill_scores()
  object <- make_gRNA_GI_for_fill(
    block_layer = "guide_pair",
    blocks = "none",
    use_blocks = TRUE,
    block_description = "old_description",
    collapse = "bio_rep"
  )

  result <- CeRberus:::fill_gRNA_GIs(object, input)

  expect_equal(result@space, object@space)
  expect_equal(result@replicates, object@replicates)
  expect_equal(result@block_layer, object@block_layer)
  expect_equal(result@blocks, object@blocks)
  expect_equal(result@use_blocks, object@use_blocks)
  expect_equal(result@block_description, object@block_description)
  expect_equal(result@collapse, object@collapse)
  expect_false(identical(result@data, object@data))
})

test_that("fill_gRNA_GIs represents missing input combinations as NA", {
  input <- make_fixed_pair_fill_scores()
  input <- input[!(
    input$gene_pair == "B;D" &
      input$guide_pair == "g2" &
      input$bio_rep == "b2" &
      input$tech_rep == "t2"
  ), ]
  object <- make_gRNA_GI_for_fill()

  result <- CeRberus:::fill_gRNA_GIs(object, input)

  expect_true(is.na(result@data["B;D", "g2", "b2", "t2"]))
  expect_equal(result@data["A;C", "g1", "b1", "t1"], 1L)
})

test_that("fill_gRNA_GIs errors for unknown value variables and missing dimensions", {
  input <- make_fixed_pair_fill_scores()
  object <- make_gRNA_GI_for_fill()

  expect_error(
    CeRberus:::fill_gRNA_GIs(object, input, value_var = "not_a_column")
  )
  expect_error(
    CeRberus:::fill_gRNA_GIs(object, input[, setdiff(names(input), "gene_pair")])
  )
  expect_error(
    CeRberus:::fill_gRNA_GIs(object, input[, setdiff(names(input), "guide_pair")])
  )
})
