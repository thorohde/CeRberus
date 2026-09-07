make_symmetric_pair_test_array <- function() {
  array(
    c(
      1, 2, 3, 4, 5, 6, 7, 8, 9,
      11, 12, 13, 14, 15, 16, 17, 18, 19
    ),
    dim = c(3L, 3L, 2L),
    dimnames = list(
      query_gene = c("A", "B", "C"),
      library_gene = c("A", "B", "C"),
      replicate = c("r1", "r2")
    )
  )
}

make_flattened_symmetric_array <- function() {
  input <- make_symmetric_pair_test_array()
  for (replicate_name in dimnames(input)[[3L]]) {
    input[, , replicate_name] <- make_symmetric(input[, , replicate_name])
  }
  input
}

test_that("flatten_symmetric_pairs returns averaged unordered pairs", {
  pairs <- c("A;A", "A;B", "A;C", "B;B", "B;C", "C;C")
  expected <- matrix(
    c(1, 11, 3, 13, 5, 15, 5, 15, 7, 17, 9, 19),
    nrow = length(pairs),
    byrow = TRUE,
    dimnames = list(pairs, c("r1", "r2"))
  )
  names(dimnames(expected)) <- c("gene_pair", "replicate")

  expect_equal(
    CeRberus:::flatten_symmetric_pairs(make_flattened_symmetric_array(), pairs),
    expected
  )
})

test_that("flatten_symmetric_pairs preserves fully missing pairs", {
  input <- make_flattened_symmetric_array()
  input["A", "B", ] <- NA_real_
  input["B", "A", ] <- NA_real_

  result <- CeRberus:::flatten_symmetric_pairs(input, pairs = "A;B")

  expected <- matrix(
    NA_real_,
    nrow = 1L,
    ncol = 2L,
    dimnames = list("A;B", c("r1", "r2"))
  )
  names(dimnames(expected)) <- c("gene_pair", "replicate")
  expect_equal(result, expected)
})

test_that("flatten_symmetric_pairs preserves one-observation shapes", {
  input <- make_flattened_symmetric_array()[, , "r1", drop = FALSE]
  cases <- list(
    all_pairs = list(
      pairs = c("A;A", "A;B", "A;C", "B;B", "B;C", "C;C"),
      expected = matrix(
        c(1, 3, 5, 5, 7, 9),
        nrow = 6L,
        ncol = 1L,
        dimnames = list(
          c("A;A", "A;B", "A;C", "B;B", "B;C", "C;C"),
          "r1"
        )
      )
    ),
    one_pair = list(
      pairs = "A;B",
      expected = matrix(3, nrow = 1L, ncol = 1L, dimnames = list("A;B", "r1"))
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expected <- case$expected
    names(dimnames(expected)) <- c("gene_pair", "replicate")
    expect_equal(
      CeRberus:::flatten_symmetric_pairs(input, case$pairs),
      expected,
      info = case_name
    )
  })
})
