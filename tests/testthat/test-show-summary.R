make_display_screen_design <- function() {
  methods::new(
    "ScreenDesign",
    query_genes = c("A", "B"),
    library_genes = c("B", "C"),
    all_genes = c("A", "B", "C"),
    query_genes_not_in_lib = "A",
    library_genes_not_in_query = "C",
    n_query_genes = 2L,
    n_lib_genes = 2L,
    n_all_genes = 3L,
    observations_per_query = c(A = 2L, B = 3L),
    all_pairs = c("A;B", "B;C"),
    unique_pairs = c("A;B", "B;C")
  )
}

make_display_guide_lfcs <- function() {
  methods::new(
    "gRNA_LFC",
    data = array(
      c(1, NA_real_, 3, Inf),
      dim = c(2L, 2L),
      dimnames = list(
        query_gene = c("A", "B"),
        library_gene = c("B", "C")
      )
    ),
    space = c("query_gene", "library_gene"),
    replicates = character(),
    collapse = "bio_rep",
    query_main_effects = c(A = 1, B = NA_real_),
    library_main_effects = c(B = 2, C = 3)
  )
}

make_display_guide_gis <- function() {
  methods::new(
    "gRNA_GI",
    data = array(
      c(1:7, NA_real_),
      dim = c(2L, 2L, 2L),
      dimnames = list(
        query_gene = c("A", "B"),
        library_gene = c("B", "C"),
        bio_rep = c("b1", "b2")
      )
    ),
    space = c("query_gene", "library_gene"),
    replicates = "bio_rep",
    block_layer = "bio_rep",
    blocks = c("1", "2"),
    use_blocks = TRUE,
    block_description = c("b1", "b2"),
    collapse = character()
  )
}

make_display_screen <- function(class = "MultiplexScreen") {
  screen <- methods::new(
    "ScreenBase",
    guideLFCs = make_display_guide_lfcs(),
    guideGIs = make_display_guide_gis(),
    screen_attr = make_display_screen_design(),
    limma_models = list(A = list(), B = list()),
    dupCorrelation = c(A = 0.1, B = 0.2),
    checks = list(
      gene_sets_equal = FALSE,
      query_sufficient = TRUE
    ),
    errors = list(
      query_genes_not_usable = "B",
      GI_computation_errors = list()
    ),
    metadata = list(
      requested_screen_type = "auto",
      inferred_screen_type = "multiplex",
      selected_screen_type = "multiplex",
      fdr_method = "BH"
    )
  )

  screen <- methods::as(screen, class)

  if (methods::is(screen, "PosAgnMultiplexScreen")) {
    screen@metadata$symmetric_analysis_method <- "preaverage"
    screen@aggregatedGuideGIs <- screen@guideGIs
    screen@aggregatedLimmaModels <- screen@limma_models
    screen@symmGeneGIs <- data.table::data.table(
      gene_pair = c("A;B", "B;C"),
      GI = c(1, -1),
      pval = c(0.01, 0.02),
      FDR = c(0.02, 0.03)
    )
  }

  screen
}

test_that("show provides compact guide-container displays", {
  lfc <- make_display_guide_lfcs()
  gi <- make_display_guide_gis()

  lfc_output <- capture.output(show(lfc))
  gi_output <- capture.output(show(gi))

  expect_match(lfc_output[[1L]], "<gRNA_LFC>", fixed = TRUE)
  expect_true(any(grepl("2 x 2", lfc_output, fixed = TRUE)))
  expect_true(any(grepl("Main effects", lfc_output, fixed = TRUE)))
  expect_false(any(grepl("@data", lfc_output, fixed = TRUE)))

  expect_match(gi_output[[1L]], "<gRNA_GI>", fixed = TRUE)
  expect_true(any(grepl("2 x 2 x 2", gi_output, fixed = TRUE)))
  expect_true(any(grepl("Blocking", gi_output, fixed = TRUE)))
  expect_false(any(grepl("@data", gi_output, fixed = TRUE)))
})

test_that("guide-container summaries return structured statistics", {
  lfc_summary <- summary(make_display_guide_lfcs())
  gi_summary <- summary(make_display_guide_gis())

  expect_s3_class(lfc_summary, "summary_gRNA_LFC")
  expect_s3_class(lfc_summary, "CeRberusSummary")
  expect_equal(lfc_summary$data$dimensions, c(2L, 2L))
  expect_equal(lfc_summary$data$values, 4L)
  expect_equal(lfc_summary$data$missing, 1L)
  expect_equal(lfc_summary$data$finite, 2L)
  expect_equal(lfc_summary$main_effects$query$total, 2L)
  expect_equal(lfc_summary$main_effects$query$missing, 1L)

  expect_s3_class(gi_summary, "summary_gRNA_GI")
  expect_s3_class(gi_summary, "CeRberusSummary")
  expect_equal(gi_summary$data$dimensions, c(2L, 2L, 2L))
  expect_equal(gi_summary$data$missing, 1L)
  expect_true(gi_summary$blocking$active)
  expect_equal(gi_summary$blocking$layer, "bio_rep")
  expect_equal(gi_summary$blocking$assignments, 2L)
})

test_that("ScreenDesign display and summary report its inferred structure", {
  design <- make_display_screen_design()
  output <- capture.output(show(design))
  result <- summary(design)

  expect_match(output[[1L]], "<ScreenDesign>", fixed = TRUE)
  expect_true(any(grepl("2 query", output, fixed = TRUE)))
  expect_true(any(grepl("2 directional", output, fixed = TRUE)))

  expect_s3_class(result, "summary_ScreenDesign")
  expect_s3_class(result, "CeRberusSummary")
  expect_equal(result$genes$query, 2L)
  expect_equal(result$genes$library, 2L)
  expect_equal(result$genes$total, 3L)
  expect_equal(result$pairs$directional, 2L)
  expect_equal(result$observations_per_query$median, 2.5)
})

test_that("screen display dispatches through the ScreenBase hierarchy", {
  classes <- c(
    "ScreenBase",
    "FixedPairScreen",
    "MultiplexScreen",
    "PosAgnMultiplexScreen"
  )

  purrr::walk(classes, function(class) {
    screen <- make_display_screen(class)
    output <- capture.output(show(screen))

    expect_match(output[[1L]], paste0("<", class, ">"), fixed = TRUE)
    expect_true(any(grepl("Genes:", output, fixed = TRUE)))
    expect_true(any(grepl("Models:", output, fixed = TRUE)))
    expect_true(any(grepl("Problems:", output, fixed = TRUE)))
    expect_false(any(grepl("@metadata", output, fixed = TRUE)))
  })
})

test_that("screen summaries reuse the screen-report schema", {
  classes <- c(
    "ScreenBase",
    "FixedPairScreen",
    "MultiplexScreen",
    "PosAgnMultiplexScreen"
  )

  purrr::walk(classes, function(class) {
    screen <- make_display_screen(class)
    result <- summary(screen)

    expect_s3_class(result, "summary_CeRberusScreen")
    expect_s3_class(result, "CeRberusSummary")
    expect_named(
      result,
      c(
        "report_version",
        "screen",
        "dimensions",
        "model",
        "checks",
        "problems",
        "results"
      )
    )
    expect_identical(result$screen$class, class)
    expect_equal(result$dimensions$query_genes, 2L)
  })
})

test_that("position-agnostic display uses symmetric results", {
  screen <- make_display_screen("PosAgnMultiplexScreen")
  output <- capture.output(show(screen))

  expect_true(any(grepl("2 symmetric gene pairs", output, fixed = TRUE)))
})

test_that("print delegates S4 objects to their show methods", {
  object <- make_display_guide_gis()

  expect_equal(capture.output(print(object)), capture.output(show(object)))
})

test_that("show and summary handle empty prototypes", {
  objects <- list(
    methods::new("gRNA_LFC"),
    methods::new("gRNA_GI"),
    methods::new("ScreenDesign"),
    methods::new("ScreenBase"),
    methods::new("FixedPairScreen"),
    methods::new("MultiplexScreen"),
    methods::new("PosAgnMultiplexScreen")
  )

  purrr::walk(objects, function(object) {
    expect_no_error(capture.output(show(object)))
    expect_no_error(summary(object))
  })
})

test_that("summary printing is readable and assignment is quiet", {
  object <- make_display_guide_gis()

  expect_output(result <- summary(object), NA)
  expect_s3_class(result, "CeRberusSummary")

  output <- capture.output(print(result))
  expect_match(output[[1L]], "<summary_gRNA_GI>", fixed = TRUE)
  expect_true(any(grepl("Data", output, fixed = TRUE)))
  expect_false(any(grepl("^List of", output)))
})
