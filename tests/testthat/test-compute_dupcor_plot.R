make_dupcor_plot_screen <- function(dupcor_data = make_dupcor_plot_data()) {
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
    metadata = list(dupcor_data = dupcor_data),
    checks = list(),
    errors = list()
  )
}

make_dupcor_plot_data <- function() {
  data.table::data.table(
    config = c("default_guide_pair_used", "default_tech_rep_used", "bio_rep_collapsed_guide_pair_used"),
    dcor = c(0.05, 0.31, -0.02),
    kept = c("selected", "", "")
  )
}

test_that("compute_dupcor_plot adds a ggplot object to every screen", {
  screens <- list(
    selected = make_dupcor_plot_screen(),
    alternative = make_dupcor_plot_screen(data.table::data.table(
      config = c("a", "b"),
      dcor = c(0.12, 0.25),
      kept = c("", "selected")
    ))
  )

  result <- CeRberus:::compute_dupcor_plot(screens)

  expect_named(result, names(screens))
  expect_s4_class(result$selected, "ScreenBase")
  expect_s3_class(result$selected@metadata$dupcor_plot, "ggplot")
  expect_s3_class(result$alternative@metadata$dupcor_plot, "ggplot")

  expect_equal(result$selected@metadata$dupcor_plot$data, screens$selected@metadata$dupcor_data)
  expect_equal(result$alternative@metadata$dupcor_plot$data, screens$alternative@metadata$dupcor_data)
  expect_equal(result$selected@metadata$dupcor_plot$labels$x, "Duplicate correlation")
  expect_equal(result$selected@metadata$dupcor_plot$labels$y, "Limma configuration")
})

test_that("compute_dupcor_plot builds the expected plot layers and fill scale", {
  result <- CeRberus:::compute_dupcor_plot(list(screen = make_dupcor_plot_screen()))
  plot <- result$screen@metadata$dupcor_plot

  expect_equal(length(plot$layers), 2L)
  expect_s3_class(plot$layers[[1L]]$geom, "GeomCol")
  expect_s3_class(plot$layers[[2L]]$geom, "GeomVline")

  vline_data <- ggplot2::ggplot_build(plot)$data[[2L]]
  expect_equal(vline_data$xintercept, c(0, 0.25))
  expect_equal(unique(vline_data$linetype), "dashed")
  expect_equal(unique(vline_data$linewidth), 1)

  fill_scale <- plot$scales$get_scales("fill")
  expect_equal(
    fill_scale$palette(2),
    setNames(c("seagreen", "grey80"), c("selected", ""))
  )
})

test_that("compute_dupcor_plot writes the plot file and creates parent directories", {
  output_file <- file.path(tempdir(), "dupcor-plot-test", "nested", "dupcor_plot.pdf")
  if (file.exists(output_file)) {
    unlink(output_file)
  }
  if (dir.exists(dirname(dirname(output_file)))) {
    unlink(dirname(dirname(output_file)), recursive = TRUE)
  }

  result <- CeRberus:::compute_dupcor_plot(
    list(screen = make_dupcor_plot_screen()),
    .fpath = output_file
  )

  expect_s3_class(result$screen@metadata$dupcor_plot, "ggplot")
  expect_true(file.exists(output_file))
  expect_gt(file.info(output_file)$size, 0)
})

test_that("compute_dupcor_plot validates input length and verbose", {
  valid_screen <- make_dupcor_plot_screen()

  expect_error(
    CeRberus:::compute_dupcor_plot(list()),
    "GI_list must contain at least one screen object"
  )
  expect_error(
    CeRberus:::compute_dupcor_plot(list(screen = valid_screen), verbose = NA),
    "verbose must be TRUE or FALSE"
  )
  expect_error(
    CeRberus:::compute_dupcor_plot(list(screen = valid_screen), verbose = c(TRUE, FALSE)),
    "verbose must be TRUE or FALSE"
  )
})
