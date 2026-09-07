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

test_that("remove_PCs supports component and preprocessing combinations", {
  x <- make_pc_matrix()
  cases <- list(
    one_pc = list(to_remove = 1, center = TRUE, scale = TRUE),
    multiple_pcs = list(to_remove = c(1, 2), center = TRUE, scale = TRUE),
    duplicate_indices = list(
      to_remove = c(1, 1, 2),
      center = TRUE,
      scale = TRUE
    ),
    no_centering = list(to_remove = 1, center = FALSE, scale = TRUE),
    no_scaling = list(to_remove = 1, center = TRUE, scale = FALSE),
    neither = list(to_remove = 1, center = FALSE, scale = FALSE)
  )

  purrr::iwalk(cases, function(case, name) {
    expect_equal(
      remove_PCs(
        x,
        to_remove = case$to_remove,
        .center = case$center,
        .scale = case$scale
      ),
      expected_remove_pcs(
        x,
        to_remove = case$to_remove,
        center = case$center,
        scale = case$scale
      ),
      tolerance = 1e-12,
      info = name
    )
  })
})

test_that("remove_PCs errors when requested PCs exceed available components", {
  x <- make_pc_matrix()

  expect_error(
    remove_PCs(x, to_remove = 4),
    "Some PCs to remove exceed the number of available principal components"
  )
})
