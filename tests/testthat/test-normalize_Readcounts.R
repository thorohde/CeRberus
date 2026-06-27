expected_normalized_readcounts <- function(readcounts, cf1 = 100, cf2 = 1) {
  total_counts <- sum(readcounts, na.rm = TRUE)
  x <- log2((readcounts / total_counts) * cf1 * length(readcounts) + cf2)
  x - min(x, na.rm = TRUE)
}

test_that("normalize_Readcounts applies library-size normalization and log transform", {
  readcounts <- c(10, 20, 40)

  result <- normalize_Readcounts(readcounts)
  expected <- expected_normalized_readcounts(readcounts)

  expect_type(result, "double")
  expect_length(result, length(readcounts))
  expect_equal(result, expected)
})

test_that("normalize_Readcounts shifts the minimum finite value to zero", {
  readcounts <- c(5, 10, 20, 40)

  result <- normalize_Readcounts(readcounts)

  expect_equal(min(result, na.rm = TRUE), 0)
  expect_true(all(result >= 0))
})

test_that("normalize_Readcounts preserves NA positions", {
  readcounts <- c(10, NA_real_, 30, 60)

  result <- normalize_Readcounts(readcounts)
  expected <- expected_normalized_readcounts(readcounts)

  expect_true(is.na(result[2]))
  expect_equal(is.na(result), is.na(readcounts))
  expect_equal(result, expected)
})

test_that("normalize_Readcounts returns all NA for zero or missing total counts", {
  expect_equal(
    normalize_Readcounts(c(0, 0, 0)),
    rep(NA_real_, 3)
  )
  expect_equal(
    normalize_Readcounts(c(NA_real_, NA_real_)),
    rep(NA_real_, 2)
  )
})

test_that("normalize_Readcounts supports custom scaling and pseudocounts", {
  readcounts <- c(1, 2, 3, 4)

  result <- normalize_Readcounts(readcounts, cf1 = 50, cf2 = 0.5)
  expected <- expected_normalized_readcounts(readcounts, cf1 = 50, cf2 = 0.5)

  expect_equal(result, expected)
})

test_that("normalize_Readcounts handles empty numeric input", {
  result <- normalize_Readcounts(numeric())

  expect_type(result, "double")
  expect_length(result, 0L)
})

test_that("normalize_Readcounts validates readcounts input", {
  expect_error(
    normalize_Readcounts(c("1", "2")),
    "readcounts must be numeric"
  )
})

test_that("normalize_Readcounts validates cf1", {
  expect_error(
    normalize_Readcounts(c(1, 2), cf1 = -1),
    "cf1 must be a single non-negative number"
  )
  expect_error(
    normalize_Readcounts(c(1, 2), cf1 = Inf),
    "cf1 must be a single non-negative number"
  )
  expect_error(
    normalize_Readcounts(c(1, 2), cf1 = c(1, 2)),
    "cf1 must be a single non-negative number"
  )
})

test_that("normalize_Readcounts validates cf2", {
  expect_error(
    normalize_Readcounts(c(1, 2), cf2 = -1),
    "cf2 must be a single non-negative number"
  )
  expect_error(
    normalize_Readcounts(c(1, 2), cf2 = NA_real_),
    "cf2 must be a single non-negative number"
  )
  expect_error(
    normalize_Readcounts(c(1, 2), cf2 = c(1, 2)),
    "cf2 must be a single non-negative number"
  )
})
