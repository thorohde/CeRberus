expected_normalized_readcounts <- function(readcounts, cf1 = 100, cf2 = 1) {
  total_counts <- sum(readcounts, na.rm = TRUE)
  x <- log2((readcounts / total_counts) * cf1 * length(readcounts) + cf2)
  x - min(x, na.rm = TRUE)
}

test_that("normalize_readcounts supports normalization options and missing values", {
  cases <- list(
    default = list(readcounts = c(10, 20, 40), cf1 = 100, cf2 = 1),
    missing = list(
      readcounts = c(10, NA_real_, 30, 60),
      cf1 = 100,
      cf2 = 1
    ),
    custom = list(readcounts = c(1, 2, 3, 4), cf1 = 50, cf2 = 0.5)
  )

  purrr::iwalk(cases, function(case, name) {
    expect_equal(
      normalize_readcounts(case$readcounts, case$cf1, case$cf2),
      expected_normalized_readcounts(
        case$readcounts,
        case$cf1,
        case$cf2
      ),
      info = name
    )
  })
})

test_that("normalize_readcounts handles degenerate inputs", {
  cases <- list(
    zero_total = list(input = c(0, 0, 0), expected = rep(NA_real_, 3)),
    all_missing = list(
      input = c(NA_real_, NA_real_),
      expected = rep(NA_real_, 2)
    ),
    empty = list(input = numeric(), expected = numeric())
  )

  purrr::iwalk(cases, function(case, name) {
    expect_equal(
      normalize_readcounts(case$input),
      case$expected,
      info = name
    )
  })
})

test_that("normalize_readcounts validates inputs and options", {
  cases <- list(
    readcounts = list(
      values = list(c("1", "2")),
      call = function(value) normalize_readcounts(value),
      error = "readcounts must be numeric"
    ),
    cf1 = list(
      values = list(-1, Inf, c(1, 2)),
      call = function(value) normalize_readcounts(c(1, 2), cf1 = value),
      error = "cf1 must be a single non-negative number"
    ),
    cf2 = list(
      values = list(-1, NA_real_, c(1, 2)),
      call = function(value) normalize_readcounts(c(1, 2), cf2 = value),
      error = "cf2 must be a single non-negative number"
    )
  )

  purrr::iwalk(cases, function(case, name) {
    purrr::walk(case$values, function(value) {
      expect_error(case$call(value), case$error, info = name)
    })
  })
})
