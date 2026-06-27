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

test_that("z_transform estimates mean and sd while ignoring missing values", {
  x <- c(1, NA_real_, 3, 5)

  result <- z_transform(x)
  expected <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

  expect_equal(result, expected)
  expect_true(is.na(result[2]))
  expect_equal(is.na(result), is.na(x))
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
