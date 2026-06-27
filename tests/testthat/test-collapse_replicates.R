make_gRNA_GI_for_collapse <- function(data, collapse = character()) {
  methods::new(
    "gRNA_GI",
    data = data,
    space = c("query_gene", "library_gene"),
    replicates = c("bio_rep", "tech_rep"),
    block_layer = character(),
    blocks = character(),
    use_blocks = FALSE,
    block_description = character(),
    collapse = collapse
  )
}

make_collapse_array <- function() {
  array(
    1:16,
    dim = c(2, 2, 2, 2),
    dimnames = list(
      query_gene = c("A", "B"),
      library_gene = c("C", "D"),
      bio_rep = c("b1", "b2"),
      tech_rep = c("t1", "t2")
    )
  )
}

test_that("collapse_replicates returns the object unchanged when no collapse layer is set", {
  data <- make_collapse_array()
  object <- make_gRNA_GI_for_collapse(data)

  result <- collapse_replicates(object)

  expect_s4_class(result, "gRNA_GI")
  expect_equal(result@data, data)
  expect_equal(result@replicates, c("bio_rep", "tech_rep"))
})

test_that("collapse_replicates collapses one replicate dimension by averaging", {
  data <- make_collapse_array()
  object <- make_gRNA_GI_for_collapse(data, collapse = "tech_rep")

  result <- collapse_replicates(object)
  expected <- apply(data, MARGIN = c(1, 2, 3), FUN = mean, na.rm = TRUE)
  expected <- array(
    expected,
    dim = c(2, 2, 2),
    dimnames = dimnames(data)[c("query_gene", "library_gene", "bio_rep")]
  )

  expect_equal(as.vector(result@data), as.vector(expected))
  expect_equal(unname(dim(result@data)), unname(dim(expected)))
  expect_equal(dimnames(result@data), dimnames(expected))
  expect_equal(result@replicates, "bio_rep")
})

test_that("collapse_replicates collapses multiple replicate dimensions", {
  data <- make_collapse_array()
  object <- make_gRNA_GI_for_collapse(
    data,
    collapse = c("bio_rep", "tech_rep")
  )

  result <- collapse_replicates(object)
  expected <- apply(data, MARGIN = c(1, 2), FUN = mean, na.rm = TRUE)
  expected <- array(
    expected,
    dim = c(2, 2),
    dimnames = dimnames(data)[c("query_gene", "library_gene")]
  )

  expect_equal(as.vector(result@data), as.vector(expected))
  expect_equal(unname(dim(result@data)), unname(dim(expected)))
  expect_equal(dimnames(result@data), dimnames(expected))
  expect_equal(result@replicates, character())
})

test_that("collapse_replicates ignores missing values when averaging", {
  data <- make_collapse_array()
  data["A", "C", "b1", "t1"] <- NA_real_
  data["A", "C", "b1", "t2"] <- 10
  object <- make_gRNA_GI_for_collapse(data, collapse = "tech_rep")

  result <- collapse_replicates(object)

  expect_equal(result@data["A", "C", "b1"], 10)
})

test_that("collapse_replicates returns NaN when all collapsed values are missing", {
  data <- make_collapse_array()
  data["A", "C", "b1", ] <- NA_real_
  object <- make_gRNA_GI_for_collapse(data, collapse = "tech_rep")

  result <- collapse_replicates(object)

  expect_true(is.nan(result@data["A", "C", "b1"]))
})

test_that("collapse_replicates preserves non-collapsed replicate dimension order", {
  data <- make_collapse_array()
  object <- make_gRNA_GI_for_collapse(data, collapse = "bio_rep")

  result <- collapse_replicates(object)

  expect_equal(
    names(dimnames(result@data)),
    c("query_gene", "library_gene", "tech_rep")
  )
  expect_equal(dimnames(result@data)$tech_rep, c("t1", "t2"))
  expect_equal(result@replicates, "tech_rep")
})

test_that("collapse_replicates rejects non-replicate collapse layers", {
  data <- make_collapse_array()
  object <- make_gRNA_GI_for_collapse(data, collapse = "query_gene")

  expect_error(
    collapse_replicates(object),
    "All collapse layers must be replicate dimensions"
  )
})

test_that("collapse_replicates reports invalid layer names", {
  data <- make_collapse_array()
  object <- make_gRNA_GI_for_collapse(
    data,
    collapse = c("tech_rep", "not_a_layer")
  )

  expect_error(
    collapse_replicates(object),
    "not_a_layer"
  )
})
