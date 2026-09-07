make_genomic_inflation_screen <- function(
  p_values,
  class = "PosAgnMultiplexScreen"
) {
  common_slots <- list(
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
    errors = list()
  )

  gene_names <- paste0("GENE", seq_along(p_values))

  if (identical(class, "FixedPairScreen")) {
    common_slots$geneGIs <- cbind(
      GI = seq_along(p_values),
      pval = p_values,
      FDR = p_values
    )
    rownames(common_slots$geneGIs) <- paste0(gene_names, ";TARGET")
  } else if (identical(class, "MultiplexScreen")) {
    common_slots$geneGIs <- array(
      data = NA_real_,
      dim = c(1L, length(p_values), 3L),
      dimnames = list(
        query_gene = "QUERY",
        library_gene = gene_names,
        variable = c("GI", "pval", "FDR")
      )
    )
    common_slots$geneGIs[,, "GI"] <- seq_along(p_values)
    common_slots$geneGIs[,, "pval"] <- p_values
    common_slots$geneGIs[,, "FDR"] <- p_values
  } else if (identical(class, "PosAgnMultiplexScreen")) {
    common_slots$symmGeneGIs <- data.table::data.table(
      gene_pair = paste0(gene_names, ";TARGET"),
      pval = p_values
    )
  } else {
    stop("Unsupported test screen class.")
  }

  do.call(methods::new, c(list(Class = class), common_slots))
}

test_that("genomic_inflation calculates lambda and ignores missing values", {
  p_values <- c(0.5, NA_real_, 0.25)
  expected <- stats::median(stats::qchisq(1 - p_values, df = 1), na.rm = TRUE) /
    stats::qchisq(0.5, df = 1)

  result <- genomic_inflation(make_genomic_inflation_screen(p_values))

  expect_equal(result, expected)
})

test_that("genomic_inflation supports result shapes and boundary values", {
  p_values <- c(0.5, 0.25)
  expected <- stats::median(stats::qchisq(1 - p_values, df = 1)) /
    stats::qchisq(0.5, df = 1)
  cases <- list(
    fixed_pair = list(class = "FixedPairScreen", p_values = p_values, expected = expected),
    multiplex = list(class = "MultiplexScreen", p_values = p_values, expected = expected),
    position_agnostic = list(
      class = "PosAgnMultiplexScreen",
      p_values = p_values,
      expected = expected
    ),
    upper_boundary = list(class = "PosAgnMultiplexScreen", p_values = 1, expected = 0),
    lower_boundary = list(class = "PosAgnMultiplexScreen", p_values = 0, expected = Inf)
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_equal(
      genomic_inflation(make_genomic_inflation_screen(case$p_values, case$class)),
      case$expected,
      info = case_name
    )
  })
})

test_that("genomic_inflation validates stored p-values", {
  missing_column <- make_genomic_inflation_screen(0.5)
  missing_column@symmGeneGIs <- data.table::data.table(gene_pair = "A;B")
  cases <- list(
    missing_column = list(
      screen = missing_column,
      error = "must contain a pval column"
    ),
    nonnumeric = list(
      screen = make_genomic_inflation_screen("0.5"),
      error = "pval column must be numeric"
    ),
    all_missing = list(
      screen = make_genomic_inflation_screen(NA_real_),
      error = "at least one non-missing p-value"
    )
  )
  for (p_value in c(-0.1, 1.1, Inf, -Inf)) {
    cases[[paste0("out_of_range_", p_value)]] <- list(
      screen = make_genomic_inflation_screen(p_value),
      error = "finite and between 0 and 1"
    )
  }

  purrr::iwalk(cases, function(case, case_name) {
    expect_error(
      genomic_inflation(case$screen),
      case$error,
      info = case_name
    )
  })
})
