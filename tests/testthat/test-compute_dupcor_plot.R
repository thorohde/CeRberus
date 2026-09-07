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
    screen_attr = methods::new("ScreenDesign"),
    dupCorrelation = numeric(),
    metadata = list(dupcor_data = dupcor_data),
    checks = list(),
    errors = list()
  )
}

make_dupcor_plot_data <- function() {
  data.table::data.table(
    config = c(
      "default_guide_pair_used",
      "default_tech_rep_used",
      "bio_rep_collapsed_guide_pair_used"
    ),
    dcor = c(0.05, 0.31, -0.02),
    kept = c("selected", "", "")
  )
}

test_that("compute_dupcor_plot adds representative plots to every screen", {
  screens <- list(
    selected = make_dupcor_plot_screen(),
    alternative = make_dupcor_plot_screen(data.table::data.table(
      config = c("a", "b"),
      dcor = c(0.12, 0.25),
      kept = c("", "selected")
    ))
  )

  result <- CeRberus:::compute_dupcor_plot(screens)
  plots <- purrr::map(result, ~ .x@metadata$dupcor_plot)
  selected_plot <- plots$selected
  vline_data <- ggplot2::ggplot_build(selected_plot)$data[[2L]]
  fill_scale <- selected_plot$scales$get_scales("fill")

  expect_equal(
    list(
      names = names(result),
      screen_classes = purrr::map_chr(result, ~ class(.x)[[1L]]),
      plot_classes = purrr::map_lgl(plots, ~ inherits(.x, "ggplot")),
      data = purrr::map(plots, "data"),
      labels = selected_plot$labels[c("x", "y")],
      geoms = unname(purrr::map_chr(
        selected_plot$layers,
        ~ class(.x$geom)[[1L]]
      )),
      thresholds = vline_data$xintercept,
      fill_palette = fill_scale$palette(2)
    ),
    list(
      names = names(screens),
      screen_classes = c(selected = "ScreenBase", alternative = "ScreenBase"),
      plot_classes = c(selected = TRUE, alternative = TRUE),
      data = purrr::map(screens, ~ .x@metadata$dupcor_data),
      labels = list(
        x = "Duplicate correlation",
        y = "Limma configuration"
      ),
      geoms = c("GeomCol", "GeomVline"),
      thresholds = c(0, 0.25),
      fill_palette = setNames(c("seagreen", "grey80"), c("selected", ""))
    )
  )
})

test_that("compute_dupcor_plot writes the plot file and creates parent directories", {
  output_file <- file.path(
    tempfile("dupcor-plot-test-"),
    "nested",
    "dupcor_plot.pdf"
  )

  result <- CeRberus:::compute_dupcor_plot(
    list(screen = make_dupcor_plot_screen()),
    .fpath = output_file
  )

  expect_equal(
    list(
      plot = inherits(result$screen@metadata$dupcor_plot, "ggplot"),
      exists = file.exists(output_file),
      nonempty = file.info(output_file)$size > 0
    ),
    list(plot = TRUE, exists = TRUE, nonempty = TRUE)
  )
})

test_that("compute_dupcor_plot validates input length and verbose", {
  valid_screen <- make_dupcor_plot_screen()
  cases <- list(
    empty = list(
      action = function() CeRberus:::compute_dupcor_plot(list()),
      error = "^gi_list must contain at least one screen object\\.$"
    ),
    missing_verbose = list(
      action = function() {
        CeRberus:::compute_dupcor_plot(
          list(screen = valid_screen),
          verbose = NA
        )
      },
      error = "^verbose must be TRUE or FALSE\\.$"
    ),
    vector_verbose = list(
      action = function() {
        CeRberus:::compute_dupcor_plot(
          list(screen = valid_screen),
          verbose = c(TRUE, FALSE)
        )
      },
      error = "^verbose must be TRUE or FALSE\\.$"
    )
  )

  purrr::iwalk(cases, function(case, name) {
    expect_error(case$action(), case$error, info = name)
  })
})
