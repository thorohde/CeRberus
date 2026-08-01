test_that("gRNA_LFC main-effect slots default to empty numeric vectors", {
  object <- methods::new("gRNA_LFC")

  expect_identical(object@query_main_effects, numeric())
  expect_identical(object@library_main_effects, numeric())
})

test_that("gRNA_LFC stores named query and library main effects", {
  object <- methods::new(
    "gRNA_LFC",
    query_main_effects = c(query_a = -0.5, query_b = 0.25),
    library_main_effects = c(library_a = 1.5, library_b = -1)
  )

  expect_identical(
    object@query_main_effects,
    c(query_a = -0.5, query_b = 0.25)
  )
  expect_identical(
    object@library_main_effects,
    c(library_a = 1.5, library_b = -1)
  )
})

test_that("compute_gene_main_effects computes positional means", {
  lfc_data <- array(
    c(1, NA, 5, NA, 3, NA, 7, NA),
    dim = c(2L, 2L, 2L),
    dimnames = list(
      query_gene = c("query_a", "query_b"),
      library_gene = c("library_a", "library_b"),
      guide_pair = c("guide_1", "guide_2")
    )
  )
  object <- methods::new(
    "gRNA_LFC",
    data = lfc_data,
    space = c("query_gene", "library_gene"),
    replicates = "guide_pair"
  )

  result <- compute_gene_main_effects(object)

  expect_identical(result@data, lfc_data)
  expect_equal(
    result@query_main_effects,
    c(query_a = 4, query_b = NA_real_)
  )
  expect_equal(
    result@library_main_effects,
    c(library_a = 2, library_b = 6)
  )
  expect_true(methods::validObject(result))
})

test_that("compute_gene_main_effects requires a directional LFC array", {
  object <- methods::new(
    "gRNA_LFC",
    data = array(
      seq_len(4L),
      dim = c(2L, 2L),
      dimnames = list(
        gene_pair = c("A;B", "A;C"),
        replicate = c("guide_1", "guide_2")
      )
    ),
    space = "gene_pair",
    replicates = "guide_pair"
  )

  expect_error(
    compute_gene_main_effects(object),
    "first two LFC dimensions must be query_gene and library_gene"
  )

  object <- methods::new(
    "gRNA_LFC",
    data = array(
      seq_len(4L),
      dim = c(2L, 2L),
      dimnames = list(
        query_gene = c("query_a", ""),
        library_gene = c("library_a", "library_b")
      )
    ),
    space = c("query_gene", "library_gene")
  )

  expect_error(
    compute_gene_main_effects(object),
    "query-gene dimension must have non-empty names"
  )
})

test_that("gRNA_LFC requires unique non-empty main-effect names", {
  expect_error(
    methods::new("gRNA_LFC", query_main_effects = c(1, 2)),
    "query_main_effects.*unique, non-empty names"
  )
  expect_error(
    methods::new(
      "gRNA_LFC",
      library_main_effects = stats::setNames(c(1, 2), c("gene", "gene"))
    ),
    "library_main_effects.*unique, non-empty names"
  )
  expect_error(
    methods::new(
      "gRNA_LFC",
      query_main_effects = stats::setNames(c(1, 2), c("gene", ""))
    ),
    "query_main_effects.*unique, non-empty names"
  )
})

test_that("gRNA_LFC requires numeric main effects", {
  expect_error(
    methods::new("gRNA_LFC", query_main_effects = c(gene = "effect")),
    "query_main_effects.*class.*numeric"
  )
})

test_that("ScreenBase validates main effects against its canonical design", {
  object <- methods::new(
    "ScreenBase",
    screen_attr = make_screen_design(
      query_genes = c("query_a", "query_b"),
      library_genes = c("library_a", "library_b")
    ),
    guideLFCs = methods::new(
      "gRNA_LFC",
      query_main_effects = c(query_b = 0.25, query_a = -0.5),
      library_main_effects = c(library_a = 1.5, library_b = -1)
    )
  )

  expect_true(methods::validObject(object))
})

test_that("ScreenBase rejects main effects inconsistent with its design", {
  expect_error(
    methods::new(
      "ScreenBase",
      screen_attr = make_screen_design(
        query_genes = c("query_a", "query_b")
      ),
      guideLFCs = methods::new(
        "gRNA_LFC",
        query_main_effects = c(query_a = -0.5)
      )
    ),
    "query_main_effects.*screen_attr\\$query_genes"
  )
  expect_error(
    methods::new(
      "ScreenBase",
      screen_attr = make_screen_design(
        library_genes = c("library_a", "library_b")
      ),
      guideLFCs = methods::new(
        "gRNA_LFC",
        library_main_effects = c(library_a = 1.5, library_c = -1)
      )
    ),
    "library_main_effects.*screen_attr\\$library_genes"
  )
})

test_that("gRNA_LFC validates array and dimension metadata", {
  expect_error(
    methods::new(
      "gRNA_LFC",
      data = array(seq_len(4L), dim = c(2L, 2L)),
      space = "query_gene"
    ),
    "array rank"
  )
  expect_error(
    methods::new(
      "gRNA_LFC",
      space = "query_gene",
      replicates = "query_gene"
    ),
    "unique dimension names"
  )
  expect_error(
    methods::new(
      "gRNA_LFC",
      space = "query_gene",
      collapse = "query_gene"
    ),
    "cannot be collapsed"
  )
})

test_that("gRNA_LFC accepts flattened replicate dimensions", {
  object <- methods::new(
    "gRNA_LFC",
    data = array(
      seq_len(4L),
      dim = c(2L, 2L),
      dimnames = list(
        query_gene = c("A", "B"),
        replicate = c("g1_b1", "g2_b1")
      )
    ),
    space = "query_gene",
    replicates = c("guide_pair", "bio_rep")
  )

  expect_true(methods::validObject(object))
})
