make_symmetric_pair_test_array <- function() {
  array(
    c(
      # replicate r1
      1,
      2,
      3,
      4,
      5,
      6,
      7,
      8,
      9,
      # replicate r2
      11,
      12,
      13,
      14,
      15,
      16,
      17,
      18,
      19
    ),
    dim = c(3L, 3L, 2L),
    dimnames = list(
      query_gene = c("A", "B", "C"),
      library_gene = c("A", "B", "C"),
      replicate = c("r1", "r2")
    )
  )
}


test_that("flatten_symmetric_pairs creates one row per unordered pair", {
  input <- make_symmetric_pair_test_array()

  for (replicate_name in dimnames(input)[[3L]]) {
    input[,, replicate_name] <- makeSymmetric(
      input[,, replicate_name]
    )
  }

  pairs <- c("A;A", "A;B", "A;C", "B;B", "B;C", "C;C")

  result <- CeRberus:::flatten_symmetric_pairs(input, pairs)

  expect_equal(dim(result), c(6L, 2L))
  expect_equal(rownames(result), pairs)
  expect_equal(colnames(result), c("r1", "r2"))
})


test_that("flatten_symmetric_pairs retains averaged orientation values", {
  original <- make_symmetric_pair_test_array()
  symmetric <- original

  for (replicate_name in dimnames(symmetric)[[3L]]) {
    symmetric[,, replicate_name] <- makeSymmetric(
      symmetric[,, replicate_name]
    )
  }

  result <- CeRberus:::flatten_symmetric_pairs(
    symmetric,
    pairs = c("A;B", "A;C", "B;C")
  )

  expect_equal(
    result["A;B", "r1"],
    mean(c(original["A", "B", "r1"], original["B", "A", "r1"]))
  )
  expect_equal(
    result["A;C", "r2"],
    mean(c(original["A", "C", "r2"], original["C", "A", "r2"]))
  )
  expect_equal(ncol(result), 2L)
})


test_that("flatten_symmetric_pairs represents self-pairs once", {
  input <- make_symmetric_pair_test_array()

  for (replicate_name in dimnames(input)[[3L]]) {
    input[,, replicate_name] <- makeSymmetric(
      input[,, replicate_name]
    )
  }

  result <- CeRberus:::flatten_symmetric_pairs(
    input,
    pairs = c("A;A", "A;B")
  )

  expect_equal(result["A;A", ], input["A", "A", ])
  expect_equal(nrow(result), 2L)
})


test_that("flatten_symmetric_pairs converts fully missing pair values to NA", {
  input <- make_symmetric_pair_test_array()
  input["A", "B", ] <- NA_real_
  input["B", "A", ] <- NA_real_

  for (replicate_name in dimnames(input)[[3L]]) {
    input[,, replicate_name] <- makeSymmetric(
      input[,, replicate_name]
    )
  }

  result <- CeRberus:::flatten_symmetric_pairs(
    input,
    pairs = "A;B"
  )

  expect_true(all(is.na(result["A;B", ])))
})


test_that("flatten_symmetric_pairs preserves matrix shape with one observation", {
  input <- make_symmetric_pair_test_array()[,, "r1", drop = FALSE]
  input[,, "r1"] <- makeSymmetric(input[,, "r1"])

  pairs <- c("A;A", "A;B", "A;C", "B;B", "B;C", "C;C")
  result <- CeRberus:::flatten_symmetric_pairs(input, pairs)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(6L, 1L))
  expect_equal(rownames(result), pairs)
  expect_equal(colnames(result), "r1")
})


test_that("flatten_symmetric_pairs handles one pair and one observation", {
  input <- make_symmetric_pair_test_array()[,, "r1", drop = FALSE]
  input[,, "r1"] <- makeSymmetric(input[,, "r1"])

  result <- CeRberus:::flatten_symmetric_pairs(input, pairs = "A;B")

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(1L, 1L))
  expect_equal(rownames(result), "A;B")
  expect_equal(colnames(result), "r1")
})
