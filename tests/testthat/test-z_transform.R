test_that("z_transform supports its standardization strategies", {
  central_x <- c(1:9, 100)
  central_bounds <- stats::quantile(central_x, probs = c(0.1, 0.9))
  central_sd <- stats::sd(central_x[
    central_x >= central_bounds[[1L]] & central_x <= central_bounds[[2L]]
  ])
  central_na_x <- c(central_x, NA_real_)
  central_na_bounds <- stats::quantile(
    central_na_x,
    probs = c(0.1, 0.9),
    na.rm = TRUE
  )
  central_na_sd <- stats::sd(
    central_na_x[
      central_na_x >= central_na_bounds[[1L]] &
        central_na_x <= central_na_bounds[[2L]]
    ],
    na.rm = TRUE
  )
  cases <- list(
    estimated = list(
      x = 1:5,
      args = list(),
      expected = (1:5 - mean(1:5)) / stats::sd(1:5)
    ),
    supplied = list(
      x = c(2, 4, 6),
      args = list(.mean = 2, .sd = 2),
      expected = c(0, 1, 2)
    ),
    central = list(
      x = central_x,
      args = list(outlier_quantile = 0.1),
      expected = (central_x - mean(central_x)) / central_sd
    ),
    central_overrides_sd = list(
      x = central_x,
      args = list(.sd = 1000, outlier_quantile = 0.1),
      expected = (central_x - mean(central_x)) / central_sd
    ),
    central_with_na = list(
      x = central_na_x,
      args = list(outlier_quantile = 0.1),
      expected = (central_na_x - mean(central_na_x, na.rm = TRUE)) /
        central_na_sd
    )
  )

  purrr::iwalk(cases, function(case, name) {
    expect_equal(
      do.call(z_transform, c(list(.x = case$x), case$args)),
      case$expected,
      info = name
    )
  })
})

test_that("z_transform preserves names and missing-value positions", {
  x <- c(a = 1, b = NA_real_, c = 3, d = 5)

  result <- z_transform(x)
  expected <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

  expect_equal(result, expected)
})

test_that("z_transform handles degenerate numeric inputs", {
  cases <- list(
    zero_sd = list(x = c(3, 3, 3), args = list(), expected = rep(NA_real_, 3)),
    zero_central_sd = list(
      x = c(0, rep(1, 8), 2),
      args = list(outlier_quantile = 0.1),
      expected = rep(NA_real_, 10)
    ),
    singleton = list(x = 42, args = list(), expected = NA_real_),
    all_missing = list(
      x = c(NA_real_, NA_real_),
      args = list(),
      expected = rep(NA_real_, 2)
    ),
    empty = list(x = numeric(), args = list(), expected = numeric())
  )

  purrr::iwalk(cases, function(case, name) {
    expect_equal(
      do.call(z_transform, c(list(.x = case$x), case$args)),
      case$expected,
      info = name
    )
  })
})

test_that("z_transform validates inputs and options", {
  cases <- list(
    input = list(
      values = list(c("1", "2"), list(1, 2), matrix(1:4, nrow = 2L)),
      call = function(value) z_transform(value),
      error = "`\\.x` must be a numeric vector"
    ),
    mean = list(
      values = list(NA_real_, Inf, c(1, 2), numeric(), "1"),
      call = function(value) z_transform(1:3, .mean = value),
      error = "`\\.mean` must be one finite numeric value"
    ),
    standard_deviation = list(
      values = list(NA_real_, Inf, -1, c(1, 2), numeric(), "1"),
      call = function(value) z_transform(1:3, .sd = value),
      error = "`\\.sd` must be one finite, non-negative numeric value"
    ),
    outlier_quantile = list(
      values = list(-0.1, 0.5, 1, Inf, NA_real_, c(0.1, 0.2), "0.1"),
      call = function(value) z_transform(1:10, outlier_quantile = value),
      error = "must be one finite numeric value in \\[0, 0.5\\)"
    )
  )

  purrr::iwalk(cases, function(case, name) {
    purrr::walk(case$values, function(value) {
      expect_error(case$call(value), case$error, info = name)
    })
  })

  expect_error(
    z_transform(1:10, .sd = -1, outlier_quantile = 0.1),
    "`\\.sd` must be one finite, non-negative numeric value"
  )
})
