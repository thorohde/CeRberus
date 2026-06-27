local_balanced_fdr <- function(pair, pval_array, method) {
  g1 <- stringr::str_split_i(pair, ";", 1)
  g2 <- stringr::str_split_i(pair, ";", 2)

  local_g1 <- c(rep(g1, ncol(pval_array)), setdiff(rownames(pval_array), g1))
  local_g2 <- c(colnames(pval_array), rep(g2, nrow(pval_array) - 1))
  local_pairs <- paste(local_g1, local_g2, sep = ";")
  local_pvals <- purrr::map2_dbl(local_g1, local_g2, ~ pval_array[.x, .y])

  stats::p.adjust(local_pvals, method = method)[local_pairs == pair]
}

test_that("balanced_FDR returns the local adjusted p-value for one pair", {
  pval_array <- matrix(
    c(
      0.90,
      0.01,
      0.20,
      0.04,
      0.80,
      0.30,
      0.50,
      0.02,
      0.70
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  result <- balanced_FDR("A;B", pval_array, fdr_method = "BH")
  expected <- local_balanced_fdr("A;B", pval_array, method = "BH")

  expect_type(result, "double")
  expect_length(result, 1L)
  expect_equal(result, expected)
})

test_that("balanced_FDR returns one value per input pair", {
  pval_array <- matrix(
    c(
      0.90,
      0.01,
      0.20,
      0.04,
      0.80,
      0.30,
      0.50,
      0.02,
      0.70
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
  pairs <- c("A;B", "B;C", "C;A")

  result <- balanced_FDR(pairs, pval_array, fdr_method = "BH")
  expected <- purrr::map_dbl(
    pairs,
    ~ local_balanced_fdr(.x, pval_array, method = "BH")
  )

  expect_length(result, length(pairs))
  expect_equal(result, expected)
})

test_that("balanced_FDR respects the requested p.adjust method", {
  pval_array <- matrix(
    c(
      0.90,
      0.01,
      0.20,
      0.04,
      0.80,
      0.30,
      0.50,
      0.02,
      0.70
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  result <- balanced_FDR("A;B", pval_array, fdr_method = "bonferroni")
  expected <- local_balanced_fdr("A;B", pval_array, method = "bonferroni")

  expect_equal(result, expected)
})

test_that("balanced_FDR treats pair direction as part of the requested result", {
  pval_array <- matrix(
    c(
      0.90,
      0.01,
      0.20,
      0.04,
      0.80,
      0.30,
      0.50,
      0.02,
      0.70
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  forward <- balanced_FDR("A;B", pval_array, fdr_method = "BH")
  reverse <- balanced_FDR("B;A", pval_array, fdr_method = "BH")

  expect_equal(forward, local_balanced_fdr("A;B", pval_array, method = "BH"))
  expect_equal(reverse, local_balanced_fdr("B;A", pval_array, method = "BH"))
  expect_false(isTRUE(all.equal(forward, reverse)))
})

test_that("balanced_FDR propagates NA p-values through p.adjust", {
  pval_array <- matrix(
    c(
      0.90,
      0.01,
      0.20,
      0.04,
      NA_real_,
      0.30,
      0.50,
      0.02,
      0.70
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  result <- balanced_FDR("B;B", pval_array, fdr_method = "BH")

  expect_true(is.na(result))
})

test_that("balanced_FDR errors for unknown genes", {
  pval_array <- matrix(
    c(0.90, 0.01, 0.04, 0.80),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("A", "B"), c("A", "B"))
  )

  expect_error(
    suppressWarnings(balanced_FDR("X;B", pval_array, fdr_method = "BH")),
    "subscript out of bounds"
  )
  expect_error(
    balanced_FDR("A;X", pval_array, fdr_method = "BH"),
    "subscript out of bounds"
  )
})

test_that("balanced_FDR errors for invalid p.adjust methods", {
  pval_array <- matrix(
    c(0.90, 0.01, 0.04, 0.80),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("A", "B"), c("A", "B"))
  )

  expect_error(
    balanced_FDR("A;B", pval_array, fdr_method = "invalid-method"),
    "match.arg"
  )
})
