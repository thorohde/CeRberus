test_that("makeSymmetric averages each matrix cell with its transpose counterpart", {
  x <- matrix(
    c(
      1,
      2,
      3,
      4,
      5,
      6,
      7,
      8,
      9
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  result <- makeSymmetric(x)
  expected <- (x + t(x)) / 2

  expect_equal(result, expected)
  expect_true(isSymmetric(result))
})

test_that("makeSymmetric preserves dimensions and dimnames", {
  x <- matrix(
    c(1, 4, 2, 5),
    nrow = 2,
    dimnames = list(c("A", "B"), c("A", "B"))
  )

  result <- makeSymmetric(x)

  expect_equal(dim(result), dim(x))
  expect_equal(rownames(result), rownames(x))
  expect_equal(colnames(result), colnames(x))
})

test_that("makeSymmetric leaves diagonal values unchanged", {
  x <- matrix(
    c(
      1,
      10,
      20,
      2
    ),
    nrow = 2,
    byrow = TRUE
  )

  result <- makeSymmetric(x)

  expect_equal(diag(result), diag(x))
})

test_that("makeSymmetric ignores one-sided missing values", {
  x <- matrix(
    c(
      1,
      NA_real_,
      3,
      4,
      5,
      NA_real_,
      NA_real_,
      8,
      9
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  result <- makeSymmetric(x)

  expect_equal(result["A", "B"], 4)
  expect_equal(result["B", "A"], 4)
  expect_equal(result["A", "C"], 3)
  expect_equal(result["C", "A"], 3)
  expect_equal(result["B", "C"], 8)
  expect_equal(result["C", "B"], 8)
})

test_that("makeSymmetric returns NaN when both directional values are missing", {
  x <- matrix(
    c(
      1,
      NA_real_,
      NA_real_,
      2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("A", "B"), c("A", "B"))
  )

  result <- makeSymmetric(x)

  expect_true(is.nan(result["A", "B"]))
  expect_true(is.nan(result["B", "A"]))
})

test_that("makeSymmetric keeps already symmetric matrices unchanged", {
  x <- matrix(
    c(
      1,
      2,
      3,
      2,
      4,
      5,
      3,
      5,
      6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  result <- makeSymmetric(x)

  expect_equal(result, x)
})
