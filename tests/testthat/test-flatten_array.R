test_that("flatten_array converts a named matrix-like array to long data.table", {
  x <- array(
    1:4,
    dim = c(2, 2),
    dimnames = list(c("A", "B"), c("C", "D"))
  )

  result <- flatten_array(x)

  expect_s3_class(result, "data.table")
  expect_named(result, c("Var1", "Var2", "value"))
  expect_equal(nrow(result), 4L)
  expect_equal(as.character(result$Var1), c("A", "B", "A", "B"))
  expect_equal(as.character(result$Var2), c("C", "C", "D", "D"))
  expect_equal(result$value, 1:4)
})

test_that("flatten_array can rename dimension columns", {
  x <- array(
    1:4,
    dim = c(2, 2),
    dimnames = list(c("A", "B"), c("C", "D"))
  )

  result <- flatten_array(x, dnames = c("query_gene", "library_gene"))

  expect_named(result, c("query_gene", "library_gene", "value"))
  expect_equal(as.character(result$query_gene), c("A", "B", "A", "B"))
  expect_equal(as.character(result$library_gene), c("C", "C", "D", "D"))
})

test_that("flatten_array supports a custom value column name", {
  x <- array(
    1:4,
    dim = c(2, 2),
    dimnames = list(c("A", "B"), c("C", "D"))
  )

  result <- flatten_array(x, value_name = "score")

  expect_named(result, c("Var1", "Var2", "score"))
  expect_equal(result$score, 1:4)
})

test_that("flatten_array generates dimension labels when dimnames are missing", {
  x <- array(1:4, dim = c(2, 2))

  result <- flatten_array(x)

  expect_named(result, c("Var1", "Var2", "value"))
  expect_equal(as.character(result$Var1), c("1", "2", "1", "2"))
  expect_equal(as.character(result$Var2), c("1", "1", "2", "2"))
  expect_equal(result$value, 1:4)
})

test_that("flatten_array supports 3D arrays", {
  x <- array(
    1:8,
    dim = c(2, 2, 2),
    dimnames = list(
      query_gene = c("A", "B"),
      library_gene = c("C", "D"),
      replicate = c("r1", "r2")
    )
  )

  result <- flatten_array(x, value_name = "score")

  expect_s3_class(result, "data.table")
  expect_named(result, c("query_gene", "library_gene", "replicate", "score"))
  expect_equal(nrow(result), 8L)
  expect_equal(result$score, c(1L, 5L, 3L, 7L, 2L, 6L, 4L, 8L))
})

test_that("flatten_array errors when input has no dim attribute", {
  expect_error(
    flatten_array(1:4),
    "Please provide an array!"
  )
})

test_that("flatten_array warns when provided dnames length does not match dimensions", {
  x <- array(
    1:4,
    dim = c(2, 2),
    dimnames = list(c("A", "B"), c("C", "D"))
  )

  expect_warning(
    result <- flatten_array(x, dnames = "query_gene"),
    "1 dimnames given, but 2 dimensions found"
  )

  expect_named(result, c("query_gene", "Var2", "value"))
})
