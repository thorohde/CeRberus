test_that("gRNA_LFC supports construction variants", {
  empty <- methods::new("gRNA_LFC")
  populated <- methods::new(
    "gRNA_LFC",
    query_main_effects = c(query_a = -0.5, query_b = 0.25),
    library_main_effects = c(library_a = 1.5, library_b = -1)
  )
  flattened <- methods::new(
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

  expect_equal(
    list(empty@query_main_effects, empty@library_main_effects),
    list(numeric(), numeric())
  )
  expect_equal(
    list(populated@query_main_effects, populated@library_main_effects),
    list(
      c(query_a = -0.5, query_b = 0.25),
      c(library_a = 1.5, library_b = -1)
    )
  )
  expect_true(methods::validObject(flattened))
})

test_that("compute_gene_main_effects computes multiplex positional means", {
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

  expect_equal(
    list(
      data = result@data,
      query = result@query_main_effects,
      library = result@library_main_effects
    ),
    list(
      data = lfc_data,
      query = c(query_a = 4, query_b = NA_real_),
      library = c(library_a = 2, library_b = 6)
    )
  )
  expect_true(methods::validObject(result))
})

test_that("compute_gene_main_effects computes fixed-pair positional means", {
  lfc_data <- array(
    c(1, 3, NA, 5, 7, NA, 9, 11, NA, 13, 15, NA),
    dim = c(3L, 2L, 2L),
    dimnames = list(
      gene_pair = c("A;C", "A;D", "B;D"),
      guide_pair = c("guide_1", "guide_2"),
      bio_rep = c("bio_1", "bio_2")
    )
  )
  object <- methods::new(
    "gRNA_LFC",
    data = lfc_data,
    space = "gene_pair",
    replicates = c("guide_pair", "bio_rep")
  )

  result <- compute_gene_main_effects(object)

  expect_equal(
    list(
      data = result@data,
      query = result@query_main_effects,
      library = result@library_main_effects
    ),
    list(
      data = lfc_data,
      query = c(A = 8, B = NA_real_),
      library = c(C = 7, D = 9)
    )
  )
  expect_true(methods::validObject(result))
})

test_that("compute_gene_main_effects validates biological dimensions", {
  cases <- list(
    malformed_pair = list(
      object = methods::new(
        "gRNA_LFC",
        data = array(
          seq_len(4L),
          dim = c(2L, 2L),
          dimnames = list(
            gene_pair = c("A;B", "invalid_pair"),
            guide_pair = c("guide_1", "guide_2")
          )
        ),
        space = "gene_pair",
        replicates = "guide_pair"
      ),
      error = "format 'query_gene;library_gene'"
    ),
    empty_query_name = list(
      object = methods::new(
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
      ),
      error = "query-gene dimension must have non-empty names"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_error(
      compute_gene_main_effects(case$object),
      case$error,
      info = case_name
    )
  })
})

test_that("gRNA_LFC rejects invalid slot and array metadata", {
  cases <- list(
    unnamed_query = list(
      call = function() methods::new("gRNA_LFC", query_main_effects = c(1, 2)),
      error = "query_main_effects.*unique, non-empty names"
    ),
    duplicate_library = list(
      call = function() methods::new(
        "gRNA_LFC",
        library_main_effects = stats::setNames(c(1, 2), c("gene", "gene"))
      ),
      error = "library_main_effects.*unique, non-empty names"
    ),
    empty_query_name = list(
      call = function() methods::new(
        "gRNA_LFC",
        query_main_effects = stats::setNames(c(1, 2), c("gene", ""))
      ),
      error = "query_main_effects.*unique, non-empty names"
    ),
    nonnumeric_effect = list(
      call = function() methods::new(
        "gRNA_LFC",
        query_main_effects = c(gene = "effect")
      ),
      error = "query_main_effects.*class.*numeric"
    ),
    wrong_array_rank = list(
      call = function() methods::new(
        "gRNA_LFC",
        data = array(seq_len(4L), dim = c(2L, 2L)),
        space = "query_gene"
      ),
      error = "array rank"
    ),
    duplicate_dimensions = list(
      call = function() methods::new(
        "gRNA_LFC",
        space = "query_gene",
        replicates = "query_gene"
      ),
      error = "unique dimension names"
    ),
    collapsed_space = list(
      call = function() methods::new(
        "gRNA_LFC",
        space = "query_gene",
        collapse = "query_gene"
      ),
      error = "cannot be collapsed"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_error(case$call(), case$error, info = case_name)
  })
})

test_that("ScreenBase validates main effects against its design", {
  valid <- methods::new(
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
  cases <- list(
    query = list(
      screen_attr = make_screen_design(query_genes = c("query_a", "query_b")),
      guide_lfcs = methods::new(
        "gRNA_LFC",
        query_main_effects = c(query_a = -0.5)
      ),
      error = "query_main_effects.*screen_attr\\$query_genes"
    ),
    library = list(
      screen_attr = make_screen_design(
        library_genes = c("library_a", "library_b")
      ),
      guide_lfcs = methods::new(
        "gRNA_LFC",
        library_main_effects = c(library_a = 1.5, library_c = -1)
      ),
      error = "library_main_effects.*screen_attr\\$library_genes"
    )
  )

  expect_true(methods::validObject(valid))
  purrr::iwalk(cases, function(case, case_name) {
    expect_error(
      methods::new(
        "ScreenBase",
        screen_attr = case$screen_attr,
        guideLFCs = case$guide_lfcs
      ),
      case$error,
      info = case_name
    )
  })
})
