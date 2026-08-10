test_that("z_transform standardizes a numeric vector with estimated mean and sd", {
  x <- c(1, 2, 3, 4, 5)

  result <- z_transform(x)
  expected <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

  expect_type(result, "double")
  expect_length(result, length(x))
  expect_equal(result, expected)
  expect_equal(mean(result), 0, tolerance = 1e-12)
  expect_equal(stats::sd(result), 1, tolerance = 1e-12)
})

test_that("z_transform supports supplied mean and standard deviation", {
  x <- c(2, 4, 6)

  result <- z_transform(x, .mean = 2, .sd = 2)

  expect_equal(result, c(0, 1, 2))
})

test_that("z_transform requires a numeric vector", {
  invalid_values <- list(
    c("1", "2"),
    list(1, 2),
    matrix(1:4, nrow = 2L)
  )

  purrr::walk(invalid_values, function(value) {
    expect_error(
      z_transform(value),
      "`\\.x` must be a numeric vector"
    )
  })
})

test_that("z_transform validates a supplied mean", {
  invalid_values <- list(
    NA_real_,
    Inf,
    c(1, 2),
    numeric(),
    "1"
  )

  purrr::walk(invalid_values, function(value) {
    expect_error(
      z_transform(1:3, .mean = value),
      "`\\.mean` must be one finite numeric value"
    )
  })
})

test_that("z_transform validates a supplied standard deviation", {
  invalid_values <- list(
    NA_real_,
    Inf,
    -1,
    c(1, 2),
    numeric(),
    "1"
  )

  purrr::walk(invalid_values, function(value) {
    expect_error(
      z_transform(1:3, .sd = value),
      "`\\.sd` must be one finite, non-negative numeric value"
    )
  })

  expect_error(
    z_transform(1:10, .sd = -1, outlier_quantile = 0.1),
    "`\\.sd` must be one finite, non-negative numeric value"
  )
})

test_that("z_transform estimates sd from a central quantile interval", {
  x <- c(1:9, 100)
  outlier_quantile <- 0.1
  bounds <- stats::quantile(
    x,
    probs = c(outlier_quantile, 1 - outlier_quantile)
  )
  central_values <- x[x >= bounds[[1L]] & x <= bounds[[2L]]]

  result <- z_transform(x, outlier_quantile = outlier_quantile)
  expected <- (x - mean(x)) / stats::sd(central_values)

  expect_equal(result, expected)
  expect_gt(result[[length(result)]], z_transform(x)[[length(x)]])
})

test_that("outlier_quantile overrides a supplied standard deviation", {
  x <- c(1:9, 100)

  result <- z_transform(x, .sd = 1000, outlier_quantile = 0.1)
  expected <- z_transform(x, outlier_quantile = 0.1)

  expect_equal(result, expected)
})

test_that("z_transform estimates mean and sd while ignoring missing values", {
  x <- c(1, NA_real_, 3, 5)

  result <- z_transform(x)
  expected <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

  expect_equal(result, expected)
  expect_true(is.na(result[2]))
  expect_equal(is.na(result), is.na(x))
})

test_that("z_transform ignores missing values when estimating central sd", {
  x <- c(1:9, 100, NA_real_)

  result <- z_transform(x, outlier_quantile = 0.1)

  expect_false(all(is.na(result[-length(result)])))
  expect_true(is.na(result[[length(result)]]))
})

test_that("z_transform preserves vector names", {
  x <- c(a = 1, b = 2, c = 3)

  result <- z_transform(x)

  expect_named(result, names(x))
})

test_that("z_transform returns all NA when standard deviation is zero", {
  result <- z_transform(c(3, 3, 3))

  expect_equal(result, rep(NA_real_, 3))
})

test_that("z_transform returns all NA when central standard deviation is zero", {
  result <- z_transform(c(0, rep(1, 8), 2), outlier_quantile = 0.1)

  expect_equal(result, rep(NA_real_, 10))
})

test_that("z_transform returns all NA when standard deviation is missing", {
  expect_equal(
    z_transform(42),
    NA_real_
  )
  expect_equal(
    z_transform(c(NA_real_, NA_real_)),
    rep(NA_real_, 2)
  )
})

test_that("z_transform handles empty numeric input", {
  result <- z_transform(numeric())

  expect_type(result, "double")
  expect_length(result, 0L)
})

test_that("z_transform validates outlier_quantile", {
  invalid_values <- list(-0.1, 0.5, 1, Inf, NA_real_, c(0.1, 0.2), "0.1")

  purrr::walk(invalid_values, function(value) {
    expect_error(
      z_transform(1:10, outlier_quantile = value),
      "must be one finite numeric value in \\[0, 0.5\\)"
    )
  })
})
