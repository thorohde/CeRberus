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

test_that("gather_symmetric_scores handles pair lookup variants", {
  scores <- make_score_matrix()
  asymmetric_scores <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("query_A", "query_B"), c("lib_X", "lib_Y", "lib_Z"))
  )
  cases <- list(
    ordinary = list(
      pairs = c("A;B", "B;C", "C;A"),
      scores = scores,
      expected = c(12, 23, 31)
    ),
    order_duplicates_and_diagonal = list(
      pairs = c("C;B", "A;A", "C;B", "B;A"),
      scores = scores,
      expected = c(32, 11, 32, 21)
    ),
    custom_separator = list(
      pairs = c("A|C", "C|B"),
      scores = scores,
      sep = "\\|",
      expected = c(13, 32)
    ),
    asymmetric_axes = list(
      pairs = c("query_A;lib_Z", "query_B;lib_X"),
      scores = asymmetric_scores,
      expected = c(3, 4)
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    args <- list(pairs = case$pairs, .arr = case$scores)
    if (!is.null(case$sep)) args$sep <- case$sep
    expect_equal(do.call(gather_symmetric_scores, args), case$expected, info = case_name)
  })
})

test_that("gather_symmetric_scores preserves missing values", {
  scores <- make_score_matrix()
  scores["A", "C"] <- NA_real_

  result <- gather_symmetric_scores(c("A;C", "C;A"), scores)

  expect_equal(result, c(NA_real_, 31))
})

test_that("gather_symmetric_scores returns an empty numeric vector for empty input", {
  scores <- make_score_matrix()

  expect_equal(gather_symmetric_scores(character(), scores), numeric())
})
