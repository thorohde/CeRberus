test_that("ScreenDesign accepts a coherent inferred design", {
  design <- make_screen_design(
    query_genes = c("A", "B"),
    library_genes = c("B", "C"),
    observations_per_query = c(A = 2L, B = 1L),
    all_pairs = c("A;B", "B;A", "B;C")
  )

  expect_true(methods::validObject(design))
  expect_equal(design@all_genes, c("A", "B", "C"))
  expect_equal(design@query_genes_not_in_lib, "A")
  expect_equal(design@library_genes_not_in_query, "C")
  expect_equal(design@unique_pairs, c("A;B", "B;C"))
})

test_that("ScreenDesign validates gene counts and derived gene sets", {
  expect_error(
    methods::new("ScreenDesign", query_genes = "A", n_query_genes = 99L),
    "n_query_genes.*length of 'query_genes'"
  )

  design <- make_screen_design(query_genes = "A", library_genes = "B")
  design@all_genes <- c("B", "A")

  expect_error(
    methods::validObject(design),
    "all_genes.*union"
  )
})

test_that("ScreenDesign validates observations per query", {
  design <- make_screen_design(
    query_genes = c("A", "B"),
    observations_per_query = c(A = 2L, B = 1L)
  )
  names(design@observations_per_query) <- c("B", "A")

  expect_error(
    methods::validObject(design),
    "observations_per_query"
  )
})

test_that("ScreenDesign validates directional and unordered pairs", {
  design <- make_screen_design(
    query_genes = c("A", "B"),
    library_genes = c("A", "B"),
    all_pairs = c("A;B", "B;A")
  )
  design@unique_pairs <- "B;A"

  expect_error(
    methods::validObject(design),
    "canonical unordered form"
  )

  design <- make_screen_design(query_genes = "A", library_genes = "B")
  design@all_pairs <- "A;missing"

  expect_error(
    methods::validObject(design),
    "genes from 'all_genes'"
  )
})
