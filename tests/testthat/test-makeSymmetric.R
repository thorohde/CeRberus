test_that("make_symmetric averages values and preserves matrix structure", {
  cases <- list(
    one_sided_missing = list(
      input = matrix(
        c(1, NA_real_, 3, 4, 5, NA_real_, NA_real_, 8, 9),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
      ),
      expected = matrix(
        c(1, 4, 3, 4, 5, 8, 3, 8, 9),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
      )
    ),
    two_sided_missing = list(
      input = matrix(
        c(1, NA_real_, NA_real_, 2),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("A", "B"), c("A", "B"))
      ),
      expected = matrix(
        c(1, NA_real_, NA_real_, 2),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("A", "B"), c("A", "B"))
      )
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_equal(make_symmetric(case$input), case$expected, info = case_name)
  })
})

test_that("make_symmetric leaves already symmetric matrices unchanged", {
  x <- matrix(
    c(1, 2, 3, 2, 4, 5, 3, 5, 6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  expect_equal(make_symmetric(x), x)
})
