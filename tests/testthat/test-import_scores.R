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
    screen_attr = list(),
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

test_that("import_scores converts valid standard input to standardized data.table", {
  input <- make_standard_import_scores()
  screen <- make_screen_for_import_scores(input)

  result <- import_scores(screen)

  expect_s4_class(result, "ScreenBase")
  expect_true(data.table::is.data.table(result@metadata$input))
  expect_true(all(c(
    "query_gene",
    "library_gene",
    "bio_rep",
    "tech_rep",
    "guide_pair",
    "GI",
    "gene_pair"
  ) %in% names(result@metadata$input)))
  expect_equal(result@metadata$input$gene_pair, c("A;D", "B;E", "C;F"))
  expect_equal(result@metadata$input$GI, c(0.1, -0.2, 0.3))
  expect_equal(result@metadata$input$extra, c("keep1", "keep2", "keep3"))
})

test_that("import_scores standardizes custom input column names", {
  input <- data.frame(
    query = c("A", "B"),
    library = c("C", "D"),
    biological = c("b1", "b2"),
    technical = c("t1", "t2"),
    guide = c("g1", "g2"),
    score = c(1.5, -2.5),
    stringsAsFactors = FALSE
  )
  screen <- make_screen_for_import_scores(
    input,
    query_col = "query",
    lib_col = "library",
    bio_rep_col = "biological",
    tech_rep_col = "technical",
    guide_col = "guide",
    gi_col = "score"
  )

  result <- import_scores(screen)

  expect_named(
    result@metadata$input,
    c("query_gene", "library_gene", "bio_rep", "tech_rep", "guide_pair", "GI", "gene_pair")
  )
  expect_equal(result@metadata$input$query_gene, c("A", "B"))
  expect_equal(result@metadata$input$library_gene, c("C", "D"))
  expect_equal(result@metadata$input$bio_rep, c("b1", "b2"))
  expect_equal(result@metadata$input$tech_rep, c("t1", "t2"))
  expect_equal(result@metadata$input$guide_pair, c("g1", "g2"))
  expect_equal(result@metadata$input$GI, c(1.5, -2.5))
  expect_equal(result@metadata$input$gene_pair, c("A;C", "B;D"))
})

test_that("import_scores copies data.table input instead of mutating it", {
  input <- data.table::as.data.table(make_standard_import_scores())
  original <- data.table::copy(input)
  screen <- make_screen_for_import_scores(input)

  result <- import_scores(screen)

  expect_equal(input, original)
  expect_false("gene_pair" %in% names(input))
  expect_true("gene_pair" %in% names(result@metadata$input))
})

test_that("import_scores tolerates absent optional replicate and score columns", {
  input <- data.frame(
    query_gene = c("A", "B"),
    library_gene = c("C", "D"),
    stringsAsFactors = FALSE
  )
  screen <- make_screen_for_import_scores(input)

  result <- import_scores(screen)

  expect_named(result@metadata$input, c("query_gene", "library_gene", "gene_pair"))
  expect_equal(result@metadata$input$gene_pair, c("A;C", "B;D"))
})

test_that("import_scores preserves metadata other than input", {
  screen <- make_screen_for_import_scores(
    make_standard_import_scores(),
    extra_metadata = list(force_fixed_pair = TRUE, custom = "keep me")
  )

  result <- import_scores(screen)

  expect_true(result@metadata$force_fixed_pair)
  expect_equal(result@metadata$custom, "keep me")
  expect_equal(result@metadata$query_col, "query_gene")
  expect_equal(result@metadata$lib_col, "library_gene")
})

test_that("import_scores validates input object and required gene columns", {
  expect_error(
    import_scores(make_screen_for_import_scores(list(query_gene = "A", library_gene = "B"))),
    "needs to be a data frame"
  )

  expect_error(
    import_scores(make_screen_for_import_scores(data.frame(library_gene = "B"))),
    "query gene column"
  )

  expect_error(
    import_scores(make_screen_for_import_scores(data.frame(query_gene = "A"))),
    "library gene column"
  )
})

test_that("import_scores validates custom required gene column names", {
  input <- data.frame(
    query = "A",
    library = "B",
    stringsAsFactors = FALSE
  )

  expect_error(
    import_scores(make_screen_for_import_scores(input, query_col = "missing_query", lib_col = "library")),
    "query gene column"
  )

  expect_error(
    import_scores(make_screen_for_import_scores(input, query_col = "query", lib_col = "missing_library")),
    "library gene column"
  )
})
