make_pval_qq_screen <- function(symmGeneGIs = make_pval_qq_data()) {
  methods::new(
    "PosAgnMultiplexScreen",
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
    metadata = list(),
    checks = list(),
    errors = list(),
    symmGeneGIs = symmGeneGIs
  )
}

make_pval_qq_data <- function() {
  data.table::data.table(
    gene_pair = c(
      "GENE1;NTC",
      "GENE2;NTC",
      "GENE1;GENE2",
      "GENE3;GENE4",
      "NTC_like;GENE5"
    ),
    pval = c(0.01, 0.20, 0.03, 0.40, 0.50)
  )
}

test_that("pvalue_qq_plot adds a ggplot object with expected labels and groups", {
  result <- CeRberus::pvalue_qq_plot(make_pval_qq_screen())

  expect_s4_class(result, "PosAgnMultiplexScreen")
  expect_s3_class(result@metadata$qq_plot, "ggplot")
  expect_equal(result@metadata$qq_plot$labels$title, "QQ-Plot")
  expect_equal(result@metadata$qq_plot$labels$x, "Expected -log10(p)")
  expect_equal(result@metadata$qq_plot$labels$y, "Observed -log10(p)")
  expect_true("qq_inflation_summary" %in% names(result@metadata))
  expect_true("qq_plot_data" %in% names(result@metadata))
  expect_setequal(
    unique(result@metadata$qq_plot$data$ctrl),
    c("target-NTC", "target-target")
  )
  expect_equal(
    result@metadata$qq_inflation_summary$n,
    c(2L, 3L)
  )
  expect_true(all(result@metadata$qq_inflation_summary$lambda >= 0))
  expect_match(result@metadata$qq_plot$labels$caption, "lambda=")
})

test_that("pvalue_qq_plot uses exact NTC matching rather than substring matching", {
  result <- CeRberus::pvalue_qq_plot(make_pval_qq_screen())

  ntc_like_group <- result@metadata$qq_plot_data[
    gene_pair == "NTC_like;GENE5",
    unique(ctrl)
  ]

  expect_identical(ntc_like_group, "target-target")
})

test_that("pvalue_qq_plot writes a file and creates parent directories", {
  output_file <- file.path(
    tempdir(),
    "qq-plot-test",
    "nested",
    "pvalue_qq_plot.pdf"
  )
  if (file.exists(output_file)) {
    unlink(output_file)
  }
  if (dir.exists(dirname(dirname(output_file)))) {
    unlink(dirname(dirname(output_file)), recursive = TRUE)
  }

  result <- CeRberus::pvalue_qq_plot(
    make_pval_qq_screen(),
    .fpath = output_file
  )

  expect_s3_class(result@metadata$qq_plot, "ggplot")
  expect_true(file.exists(output_file))
  expect_gt(file.info(output_file)$size, 0)
})

test_that("pvalue_qq_plot validates object type, verbose, and p-values", {
  invalid_screen <- methods::new(
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
    metadata = list(),
    checks = list(),
    errors = list()
  )

  expect_error(
    CeRberus::pvalue_qq_plot(invalid_screen),
    "gi_obj must be a PosAgnMultiplexScreen object"
  )
  expect_error(
    CeRberus::pvalue_qq_plot(make_pval_qq_screen(), verbose = NA),
    "verbose must be TRUE or FALSE"
  )
  expect_error(
    CeRberus::pvalue_qq_plot(
      make_pval_qq_screen(data.table::data.table(
        gene_pair = "A;B",
        pval = NA_real_
      ))
    ),
    "symmGeneGIs must contain at least one non-missing p-value"
  )
})
