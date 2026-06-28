make_multiplex_screen_for_symmetry <- function(data, replicates = dimnames(data)[[3L]]) {
  guideGIs <- methods::new(
    "gRNA_GI",
    data = data,
    space = c("query_gene", "library_gene"),
    replicates = replicates,
    block_layer = character(),
    blocks = character(),
    use_blocks = FALSE,
    block_description = character(),
    collapse = character()
  )

  methods::new(
    "MultiplexScreen",
    guideLFCs = methods::new(
      "gRNA_LFC",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character()
    ),
    guideGIs = guideGIs,
    limma_models = list(),
    geneGIs = array(numeric(), dim = 0),
    screen_attr = list(),
    dupCorrelation = numeric(),
    metadata = list(),
    checks = list(),
    errors = list()
  )
}

make_screenbase_for_symmetry <- function(data) {
  screen <- make_multiplex_screen_for_symmetry(data)
  methods::as(screen, "ScreenBase")
}

make_symmetry_array <- function(..., replicate_names = paste0("r", seq_along(list(...)))) {
  matrices <- list(...)
  stopifnot(length(matrices) > 0L)

  output <- array(
    unlist(matrices, use.names = FALSE),
    dim = c(dim(matrices[[1L]]), length(matrices)),
    dimnames = list(
      query_gene = rownames(matrices[[1L]]),
      library_gene = colnames(matrices[[1L]]),
      replicate = replicate_names
    )
  )
  output
}

make_named_matrix <- function(values) {
  matrix(
    values,
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(
      query_gene = c("G1", "G2", "G3"),
      library_gene = c("G1", "G2", "G3")
    )
  )
}

test_that("symmetry_test returns TRUE for exactly symmetric guide GI matrices", {
  symmetric_matrix <- make_named_matrix(c(
    1, 2, 3,
    2, 4, 5,
    3, 5, 6
  ))
  screen <- make_multiplex_screen_for_symmetry(
    make_symmetry_array(symmetric_matrix, replicate_names = "rep1")
  )

  expect_true(symmetry_test(screen))
})

test_that("symmetry_test accepts matrices that are symmetric within dplyr near tolerance", {
  near_symmetric_matrix <- make_named_matrix(c(
    1, 2, 3,
    2 + 1e-9, 4, 5,
    3, 5 + 1e-9, 6
  ))
  screen <- make_multiplex_screen_for_symmetry(
    make_symmetry_array(near_symmetric_matrix, replicate_names = "rep1")
  )

  expect_true(symmetry_test(screen))
})

test_that("symmetry_test uses correlation cutoff for non-exactly symmetric matrices", {
  correlated_matrix <- make_named_matrix(c(
    1.0, 2.0, 3.0,
    2.1, 4.0, 5.0,
    3.1, 5.1, 6.0
  ))
  screen <- make_multiplex_screen_for_symmetry(
    make_symmetry_array(correlated_matrix, replicate_names = "rep1")
  )

  observed_correlation <- stats::cor(
    as.vector(correlated_matrix),
    as.vector(t(correlated_matrix)),
    use = "pairwise.complete.obs"
  )

  expect_gt(observed_correlation, 0.99)
  expect_true(symmetry_test(screen, cutoff = 0.99))
  expect_false(symmetry_test(screen, cutoff = observed_correlation + 1e-6))
})

test_that("symmetry_test requires all replicate layers to pass", {
  symmetric_matrix <- make_named_matrix(c(
    1, 2, 3,
    2, 4, 5,
    3, 5, 6
  ))
  asymmetric_matrix <- make_named_matrix(c(
    1, 2, 9,
    4, 5, 6,
    7, 8, 3
  ))
  screen <- make_multiplex_screen_for_symmetry(
    make_symmetry_array(
      symmetric_matrix,
      asymmetric_matrix,
      replicate_names = c("symmetric", "asymmetric")
    )
  )

  expect_false(symmetry_test(screen, cutoff = 0.99))
})

test_that("symmetry_test uses pairwise complete observations for missing values", {
  matrix_with_missing_values <- make_named_matrix(c(
    1.0, 2.0, NA,
    2.1, 4.0, 5.0,
    NA, 5.1, 6.0
  ))
  screen <- make_multiplex_screen_for_symmetry(
    make_symmetry_array(matrix_with_missing_values, replicate_names = "rep1")
  )

  expect_true(symmetry_test(screen, cutoff = 0.99))
})

test_that("symmetry_test is only implemented for multiplex screens", {
  symmetric_matrix <- make_named_matrix(c(
    1, 2, 3,
    2, 4, 5,
    3, 5, 6
  ))
  screen <- make_screenbase_for_symmetry(
    make_symmetry_array(symmetric_matrix, replicate_names = "rep1")
  )

  expect_error(
    symmetry_test(screen),
    "unable to find an inherited method"
  )
})
