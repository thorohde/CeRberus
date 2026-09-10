make_screen_for_gi_df <- function(
  class,
  geneGIs,
  symmGeneGIs = NULL,
  guideLFCs = methods::new(
    "gRNA_LFC",
    data = array(numeric(), dim = 0),
    space = character(),
    replicates = character()
  ),
  screen_attr = methods::new("ScreenDesign")
) {
  args <- list(
    Class = class,
    guideLFCs = guideLFCs,
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
    geneGIs = geneGIs,
    screen_attr = screen_attr,
    dupCorrelation = numeric(),
    metadata = list(),
    checks = list(),
    errors = list()
  )

  if (!is.null(symmGeneGIs)) {
    args$symmGeneGIs <- symmGeneGIs
  }

  do.call(methods::new, args)
}

make_fixed_pair_geneGIs <- function() {
  matrix(
    c(
      0.1,
      0.01,
      0.02,
      -0.2,
      0.03,
      0.04
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      gene_pair = c("A;C", "B;D"),
      variable = c("GI", "pval", "FDR")
    )
  )
}

make_multiplex_geneGIs <- function() {
  array(
    c(
      0.1,
      0.2,
      -0.1,
      -0.2,
      0.01,
      0.02,
      0.03,
      0.04,
      0.05,
      0.06,
      0.07,
      0.08
    ),
    dim = c(2L, 2L, 3L),
    dimnames = list(
      c("Q1", "Q2"),
      c("L1", "L2"),
      c("GI", "pval", "FDR")
    )
  )
}

test_that("gi_df converts fixed-pair results with optional main effects", {
  geneGIs <- make_fixed_pair_geneGIs()
  guideLFCs <- methods::new(
    "gRNA_LFC",
    query_main_effects = c(A = -1.5, B = -0.5),
    library_main_effects = c(C = 0.25, D = 0.75)
  )
  screen_attr <- make_screen_design(
    query_genes = c("A", "B"),
    library_genes = c("C", "D"),
    observations_per_query = c(1L, 1L)
  )
  cases <- list(
    base = list(
      screen = make_screen_for_gi_df("FixedPairScreen", geneGIs),
      expected = data.table::data.table(
        gene_pair = c("A;C", "B;D"),
        query_gene = c("A", "B"),
        library_gene = c("C", "D"),
        GI = c(0.1, -0.2),
        pval = c(0.01, 0.03),
        FDR = c(0.02, 0.04)
      )
    ),
    main_effects = list(
      screen = make_screen_for_gi_df(
        "FixedPairScreen",
        geneGIs,
        guideLFCs = guideLFCs,
        screen_attr = screen_attr
      ),
      expected = data.table::data.table(
        gene_pair = c("A;C", "B;D"),
        query_gene = c("A", "B"),
        library_gene = c("C", "D"),
        GI = c(0.1, -0.2),
        pval = c(0.01, 0.03),
        query_main_effect = c(-1.5, -0.5),
        library_main_effect = c(0.25, 0.75),
        FDR = c(0.02, 0.04)
      )
    )
  )
  purrr::iwalk(cases, function(case, case_name) {
    expect_equal(gi_df(case$screen), case$expected, info = case_name)
  })
})

test_that("gi_df converts multiplex results with optional main effects", {
  geneGIs <- make_multiplex_geneGIs()
  guideLFCs <- methods::new(
    "gRNA_LFC",
    query_main_effects = c(Q1 = -1.5, Q2 = -0.5),
    library_main_effects = c(L1 = 0.25, L2 = 0.75)
  )
  screen_attr <- make_screen_design(
    query_genes = c("Q1", "Q2"),
    library_genes = c("L1", "L2"),
    observations_per_query = c(2L, 2L)
  )
  cases <- list(
    base = list(
      screen = make_screen_for_gi_df("MultiplexScreen", geneGIs),
      expected = data.table::data.table(
        gene_pair = c("Q1;L1", "Q1;L2", "Q2;L1", "Q2;L2"),
        query_gene = c("Q1", "Q1", "Q2", "Q2"),
        library_gene = c("L1", "L2", "L1", "L2"),
        GI = c(0.1, -0.1, 0.2, -0.2),
        pval = c(0.01, 0.03, 0.02, 0.04),
        FDR = c(0.05, 0.07, 0.06, 0.08)
      )
    ),
    main_effects = list(
      screen = make_screen_for_gi_df(
        "MultiplexScreen",
        geneGIs,
        guideLFCs = guideLFCs,
        screen_attr = screen_attr
      ),
      expected = data.table::data.table(
        gene_pair = c("Q1;L1", "Q1;L2", "Q2;L1", "Q2;L2"),
        query_gene = c("Q1", "Q1", "Q2", "Q2"),
        library_gene = c("L1", "L2", "L1", "L2"),
        GI = c(0.1, -0.1, 0.2, -0.2),
        pval = c(0.01, 0.03, 0.02, 0.04),
        query_main_effect = c(-1.5, -1.5, -0.5, -0.5),
        library_main_effect = c(0.25, 0.75, 0.25, 0.75),
        FDR = c(0.05, 0.07, 0.06, 0.08)
      )
    )
  )
  data.table::setkey(cases$base$expected, query_gene, library_gene)
  data.table::setkey(cases$main_effects$expected, query_gene, library_gene)

  purrr::iwalk(cases, function(case, case_name) {
    expect_equal(gi_df(case$screen), case$expected, info = case_name)
  })
})

test_that("gi_df returns symmetrized data for position-agnostic multiplex screens", {
  symmGeneGIs <- data.table::data.table(
    gene_pair = c("A;B", "A;C"),
    query_gene = c("A", "A"),
    library_gene = c("B", "C"),
    GI = c(0.4, -0.5),
    GI_z = c(1.2, -1.4),
    pval = c(0.01, 0.02),
    FDR = c(0.03, 0.04),
    GI_ab = c(0.3, -0.4),
    pval_ab = c(0.02, 0.03),
    FDR_ab = c(0.04, 0.05),
    GI_ba = c(0.5, -0.6),
    pval_ba = c(0.04, 0.05),
    FDR_ba = c(0.06, 0.07),
    GI_aggregated = c(0.4, -0.5),
    pval_aggregated = c(0.01, 0.02),
    FDR_aggregated = c(0.03, 0.04)
  )
  screen <- make_screen_for_gi_df(
    "PosAgnMultiplexScreen",
    make_multiplex_geneGIs(),
    symmGeneGIs = symmGeneGIs
  )

  result <- gi_df(screen)

  expect_s3_class(result, "data.table")
  expect_equal(result, symmGeneGIs)
  expect_equal(result$GI, result$GI_aggregated)
  expect_equal(result$pval, result$pval_aggregated)
  expect_equal(result$FDR, result$FDR_aggregated)
})

test_that("gi_df appends stored NTC annotations for every result shape", {
  symmetric <- data.table::data.table(
    gene_pair = c("A;B", "A;C"),
    query_gene = c("A", "A"),
    library_gene = c("B", "C")
  )
  cases <- list(
    fixed = make_screen_for_gi_df("FixedPairScreen", make_fixed_pair_geneGIs()),
    multiplex = make_screen_for_gi_df("MultiplexScreen", make_multiplex_geneGIs()),
    symmetric = make_screen_for_gi_df(
      "PosAgnMultiplexScreen",
      make_multiplex_geneGIs(),
      symmGeneGIs = symmetric
    )
  )

  purrr::iwalk(cases, function(screen, shape) {
    pairs <- gi_df(screen)$gene_pair
    screen@metadata$ntc_pair_annotations <- data.table::data.table(
      gene_pair = pairs,
      has_NTC = seq_along(pairs) %% 2L == 1L
    )

    expect_identical(
      gi_df(screen)$has_NTC,
      seq_along(pairs) %% 2L == 1L,
      info = shape
    )
  })
})

test_that("gi_df preserves row order from fixed-pair geneGIs row names", {
  geneGIs <- make_fixed_pair_geneGIs()
  geneGIs <- geneGIs[c("B;D", "A;C"), ]
  screen <- make_screen_for_gi_df("FixedPairScreen", geneGIs)

  result <- gi_df(screen)

  expect_equal(result$gene_pair, c("B;D", "A;C"))
  expect_equal(result$query_gene, c("B", "A"))
  expect_equal(result$library_gene, c("D", "C"))
})

test_that("gi_df leaves missing fixed-pair library genes as NA when row names lack separators", {
  geneGIs <- make_fixed_pair_geneGIs()
  rownames(geneGIs) <- c("AC", "BD")
  screen <- make_screen_for_gi_df("FixedPairScreen", geneGIs)

  result <- gi_df(screen)

  expect_equal(result$query_gene, c("AC", "BD"))
  expect_true(all(is.na(result$library_gene)))
})
