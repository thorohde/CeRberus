test_that("ScreenDesign accepts and derives a coherent design", {
  design <- make_screen_design(
    query_genes = c("A", "B"),
    library_genes = c("B", "C"),
    observations_per_query = c(A = 2L, B = 1L),
    all_pairs = c("A;B", "B;A", "B;C")
  )

  expect_true(methods::validObject(design))
  expect_equal(
    list(
      all_genes = design@all_genes,
      query_missing = design@query_genes_not_in_lib,
      library_missing = design@library_genes_not_in_query,
      unique_pairs = design@unique_pairs
    ),
    list(
      all_genes = c("A", "B", "C"),
      query_missing = "A",
      library_missing = "C",
      unique_pairs = c("A;B", "B;C")
    )
  )
})

test_that("ScreenDesign rejects invalid slot mutations", {
  cases <- list(
    gene_count = list(
      call = function() methods::new(
        "ScreenDesign",
        query_genes = "A",
        n_query_genes = 99L
      ),
      error = "n_query_genes.*length of 'query_genes'"
    ),
    all_genes = list(
      call = function() {
        design <- make_screen_design(query_genes = "A", library_genes = "B")
        design@all_genes <- c("B", "A")
        methods::validObject(design)
      },
      error = "all_genes.*union"
    ),
    observations = list(
      call = function() {
        design <- make_screen_design(
          query_genes = c("A", "B"),
          observations_per_query = c(A = 2L, B = 1L)
        )
        names(design@observations_per_query) <- c("B", "A")
        methods::validObject(design)
      },
      error = "observations_per_query"
    ),
    unordered_pairs = list(
      call = function() {
        design <- make_screen_design(
          query_genes = c("A", "B"),
          library_genes = c("A", "B"),
          all_pairs = c("A;B", "B;A")
        )
        design@unique_pairs <- "B;A"
        methods::validObject(design)
      },
      error = "canonical unordered form"
    ),
    unknown_pair_gene = list(
      call = function() {
        design <- make_screen_design(query_genes = "A", library_genes = "B")
        design@all_pairs <- "A;missing"
        methods::validObject(design)
      },
      error = "genes from 'all_genes'"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_error(case$call(), case$error, info = case_name)
  })
})
