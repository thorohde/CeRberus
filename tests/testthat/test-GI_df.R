make_screen_for_gi_df <- function(class, geneGIs, symmGeneGIs = NULL) {
  args <- list(
    Class = class,
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
    geneGIs = geneGIs,
    screen_attr = list(),
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

test_that("gi_df converts fixed-pair gene GI matrices to data.table output", {
  geneGIs <- make_fixed_pair_geneGIs()
  screen <- make_screen_for_gi_df("FixedPairScreen", geneGIs)

  result <- gi_df(screen)

  expect_s3_class(result, "data.table")
  expect_named(
    result,
    c("gene_pair", "query_gene", "library_gene", "GI", "pval", "FDR")
  )
  expect_equal(result$gene_pair, c("A;C", "B;D"))
  expect_equal(result$query_gene, c("A", "B"))
  expect_equal(result$library_gene, c("C", "D"))
  expect_equal(result$GI, c(0.1, -0.2))
  expect_equal(result$pval, c(0.01, 0.03))
  expect_equal(result$FDR, c(0.02, 0.04))
})

test_that("gi_df converts multiplex gene GI arrays to long-wide data.table output", {
  geneGIs <- make_multiplex_geneGIs()
  screen <- make_screen_for_gi_df("MultiplexScreen", geneGIs)

  result <- gi_df(screen)

  expect_s3_class(result, "data.table")
  expect_named(
    result,
    c("gene_pair", "query_gene", "library_gene", "GI", "pval", "FDR")
  )
  expect_equal(nrow(result), 4L)
  expect_equal(
    result$gene_pair,
    paste(result$query_gene, result$library_gene, sep = ";")
  )
  expect_equal(result$query_gene, c("Q1", "Q1", "Q2", "Q2"))
  expect_equal(result$library_gene, c("L1", "L2", "L1", "L2"))

  expected <- data.table::data.table(
    query_gene = c("Q1", "Q1", "Q2", "Q2"),
    library_gene = c("L1", "L2", "L1", "L2"),
    GI = c(
      geneGIs["Q1", "L1", "GI"],
      geneGIs["Q1", "L2", "GI"],
      geneGIs["Q2", "L1", "GI"],
      geneGIs["Q2", "L2", "GI"]
    ),
    pval = c(
      geneGIs["Q1", "L1", "pval"],
      geneGIs["Q1", "L2", "pval"],
      geneGIs["Q2", "L1", "pval"],
      geneGIs["Q2", "L2", "pval"]
    ),
    FDR = c(
      geneGIs["Q1", "L1", "FDR"],
      geneGIs["Q1", "L2", "FDR"],
      geneGIs["Q2", "L1", "FDR"],
      geneGIs["Q2", "L2", "FDR"]
    )
  )
  expected[, gene_pair := paste(query_gene, library_gene, sep = ";")]
  expected <- expected[, .SD, .SDcols = names(result)]

  expect_equal(as.data.frame(result), as.data.frame(expected))
})

test_that("gi_df returns symmetrized data for position-agnostic multiplex screens", {
  symmGeneGIs <- data.table::data.table(
    gene_pair = c("A;B", "A;C"),
    query_gene = c("A", "A"),
    library_gene = c("B", "C"),
    GI = c(0.4, -0.5),
    GI_z = c(1.2, -1.4),
    pval = c(0.01, 0.02),
    FDR = c(0.03, 0.04)
  )
  screen <- make_screen_for_gi_df(
    "PosAgnMultiplexScreen",
    make_multiplex_geneGIs(),
    symmGeneGIs = symmGeneGIs
  )

  result <- gi_df(screen)

  expect_s3_class(result, "data.table")
  expect_equal(result, symmGeneGIs)
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
