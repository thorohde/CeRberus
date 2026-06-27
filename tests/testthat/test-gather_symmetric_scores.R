make_score_matrix <- function() {
  matrix(
    c(
      11, 12, 13,
      21, 22, 23,
      31, 32, 33
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
}

test_that("gather_symmetric_scores returns matrix values for requested pairs", {
  scores <- make_score_matrix()

  result <- gather_symmetric_scores(c("A;B", "B;C", "C;A"), scores)

  expect_type(result, "double")
  expect_equal(result, c(12, 23, 31))
})

test_that("gather_symmetric_scores preserves input order and duplicates", {
  scores <- make_score_matrix()
  pairs <- c("C;B", "A;A", "C;B", "B;A")

  result <- gather_symmetric_scores(pairs, scores)

  expect_length(result, length(pairs))
  expect_equal(result, c(32, 11, 32, 21))
})

test_that("gather_symmetric_scores treats pair direction as a matrix lookup", {
  scores <- make_score_matrix()

  forward <- gather_symmetric_scores("A;B", scores)
  reverse <- gather_symmetric_scores("B;A", scores)

  expect_equal(forward, scores["A", "B"])
  expect_equal(reverse, scores["B", "A"])
  expect_false(isTRUE(all.equal(forward, reverse)))
})

test_that("gather_symmetric_scores returns diagonal self-pair scores", {
  scores <- make_score_matrix()

  result <- gather_symmetric_scores(c("A;A", "B;B", "C;C"), scores)

  expect_equal(result, c(11, 22, 33))
})

test_that("gather_symmetric_scores supports custom pair separators", {
  scores <- make_score_matrix()

  result <- gather_symmetric_scores(c("A|C", "C|B"), scores, sep = "\\|")

  expect_equal(result, c(13, 32))
})

test_that("gather_symmetric_scores preserves missing values from the score matrix", {
  scores <- make_score_matrix()
  scores["A", "C"] <- NA_real_

  result <- gather_symmetric_scores(c("A;C", "C;A"), scores)

  expect_true(is.na(result[1]))
  expect_equal(result[2], 31)
})

test_that("gather_symmetric_scores works with asymmetric row and column gene sets", {
  scores <- matrix(
    c(
      1, 2, 3,
      4, 5, 6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("query_A", "query_B"), c("lib_X", "lib_Y", "lib_Z"))
  )

  result <- gather_symmetric_scores(
    c("query_A;lib_Z", "query_B;lib_X"),
    scores
  )

  expect_equal(result, c(3, 4))
})

test_that("gather_symmetric_scores validates that first genes are row names", {
  scores <- make_score_matrix()

  expect_error(
    gather_symmetric_scores("X;A", scores),
    "length\\(setdiff\\(genes1, rownames\\(.arr\\)\\)\\) == 0 is not TRUE"
  )
})

test_that("gather_symmetric_scores validates that second genes are column names", {
  scores <- make_score_matrix()

  expect_error(
    gather_symmetric_scores("A;X", scores),
    "length\\(setdiff\\(genes2, colnames\\(.arr\\)\\)\\) == 0 is not TRUE"
  )
})

test_that("gather_symmetric_scores returns an empty numeric vector for empty input", {
  scores <- make_score_matrix()

  result <- gather_symmetric_scores(character(), scores)

  expect_type(result, "double")
  expect_length(result, 0L)
})
