make_pc_matrix <- function() {
  matrix(
    c(
      1,
      2,
      4,
      2,
      3,
      7,
      3,
      5,
      9,
      4,
      7,
      12,
      5,
      11,
      15
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      paste0("sample", 1:5),
      paste0("feature", 1:3)
    )
  )
}

expected_remove_pcs <- function(x, to_remove, center = TRUE, scale = TRUE) {
  pca_result <- stats::prcomp(x, center = center, scale. = scale)
  pcs <- pca_result$x
  pcs[, unique(to_remove)] <- 0

  reconstructed <- pcs %*% t(pca_result$rotation)
  if (scale) {
    reconstructed <- sweep(reconstructed, 2, pca_result$scale, FUN = "*")
  }
  if (center) {
    reconstructed <- sweep(reconstructed, 2, pca_result$center, FUN = "+")
  }
  reconstructed
}

test_that("remove_PCs returns input unchanged when to_remove is NA", {
  x <- make_pc_matrix()

  result <- remove_PCs(x, to_remove = NA)

  expect_equal(result, x)
})

test_that("remove_PCs removes a selected principal component", {
  x <- make_pc_matrix()

  result <- remove_PCs(x, to_remove = 1)
  expected <- expected_remove_pcs(x, to_remove = 1)

  expect_type(result, "double")
  expect_equal(unname(dim(result)), unname(dim(x)))
  expect_equal(result, expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(result, x)))
})

test_that("remove_PCs can remove multiple principal components", {
  x <- make_pc_matrix()

  result <- remove_PCs(x, to_remove = c(1, 2))
  expected <- expected_remove_pcs(x, to_remove = c(1, 2))

  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("remove_PCs ignores duplicate PC indices", {
  x <- make_pc_matrix()

  result_with_duplicates <- remove_PCs(x, to_remove = c(1, 1, 2))
  result_without_duplicates <- remove_PCs(x, to_remove = c(1, 2))

  expect_equal(result_with_duplicates, result_without_duplicates)
})

test_that("remove_PCs supports disabling centering", {
  x <- make_pc_matrix()

  result <- remove_PCs(x, to_remove = 1, .center = FALSE, .scale = TRUE)
  expected <- expected_remove_pcs(
    x,
    to_remove = 1,
    center = FALSE,
    scale = TRUE
  )

  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("remove_PCs supports disabling scaling", {
  x <- make_pc_matrix()

  result <- remove_PCs(x, to_remove = 1, .center = TRUE, .scale = FALSE)
  expected <- expected_remove_pcs(
    x,
    to_remove = 1,
    center = TRUE,
    scale = FALSE
  )

  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("remove_PCs can remove PCs without centering or scaling", {
  x <- make_pc_matrix()

  result <- remove_PCs(x, to_remove = 1, .center = FALSE, .scale = FALSE)
  expected <- expected_remove_pcs(
    x,
    to_remove = 1,
    center = FALSE,
    scale = FALSE
  )

  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("remove_PCs errors when requested PCs exceed available components", {
  x <- make_pc_matrix()

  expect_error(
    remove_PCs(x, to_remove = 4),
    "Some PCs to remove exceed the number of available principal components"
  )
})

test_that("remove_PCs propagates prcomp errors for invalid input", {
  expect_error(
    remove_PCs(matrix(c("a", "b", "c", "d"), nrow = 2), to_remove = 1)
  )
})
