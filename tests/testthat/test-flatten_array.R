test_that("flatten_array supports 2D naming options", {
  x <- array(
    1:4,
    dim = c(2, 2),
    dimnames = list(c("A", "B"), c("C", "D"))
  )
  expected <- data.table::data.table(
    Var1 = factor(c("A", "B", "A", "B"), levels = c("A", "B")),
    Var2 = factor(c("C", "C", "D", "D"), levels = c("C", "D")),
    value = 1:4
  )
  renamed <- data.table::copy(expected)
  data.table::setnames(renamed, c("Var1", "Var2"), c("query_gene", "library_gene"))
  custom_value <- data.table::copy(expected)
  data.table::setnames(custom_value, "value", "score")
  generated_labels <- data.table::data.table(
    Var1 = factor(c("1", "2", "1", "2"), levels = c("1", "2")),
    Var2 = factor(c("1", "1", "2", "2"), levels = c("1", "2")),
    value = 1:4
  )
  cases <- list(
    default = list(action = function() flatten_array(x), expected = expected),
    renamed = list(
      action = function() {
        flatten_array(x, dnames = c("query_gene", "library_gene"))
      },
      expected = renamed
    ),
    custom_value = list(
      action = function() flatten_array(x, value_name = "score"),
      expected = custom_value
    ),
    generated_labels = list(
      action = function() flatten_array(array(1:4, dim = c(2, 2))),
      expected = generated_labels
    )
  )

  purrr::iwalk(cases, function(case, name) {
    expect_equal(case$action(), case$expected, info = name)
  })
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
  expected <- data.table::data.table(
    query_gene = rep(c("A", "B"), each = 4),
    library_gene = rep(rep(c("C", "D"), each = 2), 2),
    replicate = rep(c("r1", "r2"), 4),
    score = c(1L, 5L, 3L, 7L, 2L, 6L, 4L, 8L)
  )
  data.table::setkey(expected, query_gene, library_gene, replicate)

  expect_equal(result, expected)
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

  expect_equal(
    result,
    data.table::data.table(
      query_gene = factor(c("A", "B", "A", "B"), levels = c("A", "B")),
      Var2 = factor(c("C", "C", "D", "D"), levels = c("C", "D")),
      value = 1:4
    )
  )
})
