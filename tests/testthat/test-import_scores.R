make_screen_for_import_scores <- function(
  input,
  query_col = "query_gene",
  lib_col = "library_gene",
  bio_rep_col = "bio_rep",
  tech_rep_col = "tech_rep",
  guide_col = "guide_pair",
  gi_col = "GI",
  extra_metadata = list()
) {
  metadata <- c(
    list(
      input = input,
      query_col = query_col,
      lib_col = lib_col,
      bio_rep_col = bio_rep_col,
      tech_rep_col = tech_rep_col,
      guide_col = guide_col,
      gi_col = gi_col
    ),
    extra_metadata
  )

  methods::new(
    "ScreenBase",
    guideLFCs = methods::new(
      "gRNA_LFC",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character()
    ),
    guideGIs = methods::new(
      "gRNA_GI",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character(),
      block_layer = character(),
      blocks = character(),
      use_blocks = FALSE,
      block_description = character(),
      collapse = character()
    ),
    limma_models = list(),
    geneGIs = array(numeric(), dim = 0),
    screen_attr = methods::new("ScreenDesign"),
    dupCorrelation = numeric(),
    metadata = metadata,
    checks = list(),
    errors = list()
  )
}

make_standard_import_scores <- function() {
  data.frame(
    query_gene = c("A", "B", "C"),
    library_gene = c("D", "E", "F"),
    bio_rep = c("b1", "b1", "b2"),
    tech_rep = c("t1", "t2", "t1"),
    guide_pair = c("g1", "g2", "g3"),
    GI = c(0.1, -0.2, 0.3),
    extra = c("keep1", "keep2", "keep3"),
    stringsAsFactors = FALSE
  )
}

test_that("import_scores standardizes supported input variants", {
  custom_input <- data.frame(
    query = c("A", "B"),
    library = c("C", "D"),
    biological = c("b1", "b2"),
    technical = c("t1", "t2"),
    guide = c("g1", "g2"),
    score = c(1.5, -2.5),
    stringsAsFactors = FALSE
  )
  minimal_input <- data.frame(
    query_gene = c("A", "B"),
    library_gene = c("C", "D"),
    stringsAsFactors = FALSE
  )
  cases <- list(
    standard = list(
      screen = make_screen_for_import_scores(make_standard_import_scores()),
      expected = data.table::data.table(
        query_gene = c("A", "B", "C"),
        library_gene = c("D", "E", "F"),
        bio_rep = c("b1", "b1", "b2"),
        tech_rep = c("t1", "t2", "t1"),
        guide_pair = c("g1", "g2", "g3"),
        GI = c(0.1, -0.2, 0.3),
        extra = c("keep1", "keep2", "keep3"),
        gene_pair = c("A;D", "B;E", "C;F")
      )
    ),
    custom_columns = list(
      screen = make_screen_for_import_scores(
        custom_input,
        query_col = "query",
        lib_col = "library",
        bio_rep_col = "biological",
        tech_rep_col = "technical",
        guide_col = "guide",
        gi_col = "score"
      ),
      expected = data.table::data.table(
        query_gene = c("A", "B"),
        library_gene = c("C", "D"),
        bio_rep = c("b1", "b2"),
        tech_rep = c("t1", "t2"),
        guide_pair = c("g1", "g2"),
        GI = c(1.5, -2.5),
        gene_pair = c("A;C", "B;D")
      )
    ),
    optional_columns_absent = list(
      screen = make_screen_for_import_scores(minimal_input),
      expected = data.table::data.table(
        query_gene = c("A", "B"),
        library_gene = c("C", "D"),
        gene_pair = c("A;C", "B;D")
      )
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    result <- import_scores(case$screen)

    expect_equal(result@metadata$input, case$expected, info = case_name)
  })
})

test_that("import_scores preserves caller input and unrelated metadata", {
  input <- data.table::as.data.table(make_standard_import_scores())
  original <- data.table::copy(input)
  screen <- make_screen_for_import_scores(
    input,
    extra_metadata = list(force_fixed_pair = TRUE, custom = "keep me")
  )

  result <- import_scores(screen)

  expect_equal(input, original)
  expect_false("gene_pair" %in% names(input))
  expect_true("gene_pair" %in% names(result@metadata$input))
  expect_equal(
    result@metadata[c("force_fixed_pair", "custom", "query_col", "lib_col")],
    list(
      force_fixed_pair = TRUE,
      custom = "keep me",
      query_col = "query_gene",
      lib_col = "library_gene"
    )
  )
})

test_that("import_scores validates input and required gene columns", {
  input <- data.frame(
    query = "A",
    library = "B",
    stringsAsFactors = FALSE
  )
  cases <- list(
    non_data_frame = list(
      screen = make_screen_for_import_scores(list(
        query_gene = "A",
        library_gene = "B"
      )),
      error = "needs to be a data frame"
    ),
    missing_query = list(
      screen = make_screen_for_import_scores(data.frame(library_gene = "B")),
      error = "query gene column"
    ),
    missing_library = list(
      screen = make_screen_for_import_scores(data.frame(query_gene = "A")),
      error = "library gene column"
    ),
    missing_custom_query = list(
      screen = make_screen_for_import_scores(
        input,
        query_col = "missing_query",
        lib_col = "library"
      ),
      error = "query gene column"
    ),
    missing_custom_library = list(
      screen = make_screen_for_import_scores(
        input,
        query_col = "query",
        lib_col = "missing_library"
      ),
      error = "library gene column"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_error(import_scores(case$screen), case$error, info = case_name)
  })
})

test_that("import_scores rejects unsupported columns", {
  input <- make_standard_import_scores()
  cases <- list(
    condition = list(
      input = transform(input, condition = "treated"),
      error = "Unsupported input column\\(s\\): condition"
    ),
    orientation = list(
      input = transform(input, orientation = "forward"),
      error = "Unsupported input column\\(s\\): orientation"
    ),
    both = list(
      input = transform(input, condition = "treated", orientation = "forward"),
      error = "Unsupported input column\\(s\\): condition, orientation"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_error(
      import_scores(make_screen_for_import_scores(case$input)),
      case$error,
      info = case_name
    )
  })
})
