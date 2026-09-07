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
    screen_attr = methods::new("ScreenDesign"),
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

test_that("pvalue_qq_plot stores representative plot data and diagnostics", {
  result <- CeRberus::pvalue_qq_plot(make_pval_qq_screen())
  plot <- result@metadata$qq_plot
  plot_data <- result@metadata$qq_plot_data
  summary <- result@metadata$qq_inflation_summary

  expect_equal(
    list(
      screen_class = class(result)[[1L]],
      plot_class = inherits(plot, "ggplot"),
      labels = plot$labels[c("title", "x", "y")],
      metadata = c(
        qq_plot_data = "qq_plot_data" %in% names(result@metadata),
        qq_inflation_summary = "qq_inflation_summary" %in%
          names(result@metadata)
      ),
      groups = sort(unique(plot_data$ctrl)),
      ntc_like_group = plot_data[
        gene_pair == "NTC_like;GENE5",
        unique(ctrl)
      ],
      group_sizes = stats::setNames(summary$n, summary$ctrl),
      nonnegative_inflation = all(summary$lambda >= 0),
      caption_has_inflation = grepl("lambda=", plot$labels$caption)
    ),
    list(
      screen_class = "PosAgnMultiplexScreen",
      plot_class = TRUE,
      labels = list(
        title = "QQ-Plot",
        x = "Expected -log10(p)",
        y = "Observed -log10(p)"
      ),
      metadata = c(qq_plot_data = TRUE, qq_inflation_summary = TRUE),
      groups = c("target-NTC", "target-target"),
      ntc_like_group = "target-target",
      group_sizes = c(`target-NTC` = 2L, `target-target` = 3L),
      nonnegative_inflation = TRUE,
      caption_has_inflation = TRUE
    )
  )
})

test_that("pvalue_qq_plot writes a file and creates parent directories", {
  output_file <- file.path(
    tempfile("qq-plot-test-"),
    "nested",
    "pvalue_qq_plot.pdf"
  )

  result <- CeRberus::pvalue_qq_plot(
    make_pval_qq_screen(),
    .fpath = output_file
  )

  expect_equal(
    list(
      plot = inherits(result@metadata$qq_plot, "ggplot"),
      exists = file.exists(output_file),
      nonempty = file.info(output_file)$size > 0
    ),
    list(plot = TRUE, exists = TRUE, nonempty = TRUE)
  )
})

test_that("pvalue_qq_plot validates object type, verbose, and p-values", {
  cases <- list(
    object = list(
      action = function() CeRberus::pvalue_qq_plot(methods::new("ScreenBase")),
      error = "^gi_obj must be a PosAgnMultiplexScreen object\\.$"
    ),
    verbose = list(
      action = function() {
        CeRberus::pvalue_qq_plot(make_pval_qq_screen(), verbose = NA)
      },
      error = "^verbose must be TRUE or FALSE\\.$"
    ),
    p_values = list(
      action = function() {
        CeRberus::pvalue_qq_plot(
          make_pval_qq_screen(data.table::data.table(
            gene_pair = "A;B",
            pval = NA_real_
          ))
        )
      },
      error = paste0(
        "^symmGeneGIs must contain at least one non-missing p-value\\.$"
      )
    )
  )

  purrr::iwalk(cases, function(case, name) {
    expect_error(case$action(), case$error, info = name)
  })
})
