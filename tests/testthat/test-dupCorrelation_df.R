make_screen_for_dup_correlation_df <- function(dupcor) {
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
    dupCorrelation = dupcor,
    metadata = list(),
    checks = list(),
    errors = list()
  )
}

test_that("dup_correlation_df converts correlation vectors without modifying input", {
  cases <- list(
    scalar = 0.123,
    named_vector = c(Q1 = 0.1, Q2 = 0.2, Q3 = NA_real_),
    empty = numeric()
  )

  purrr::iwalk(cases, function(dupcor, name) {
    screen <- make_screen_for_dup_correlation_df(dupcor)
    original <- screen

    result <- dup_correlation_df(screen)

    expect_equal(
      list(result = result, input = screen),
      list(
        result = data.table::data.table(dupcor = unname(dupcor)),
        input = original
      ),
      info = name
    )
  })
})
