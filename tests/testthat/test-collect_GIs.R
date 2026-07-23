make_collect_model <- function(genes, coefficients, pvalues) {
  list(
    coefficients = matrix(
      coefficients,
      ncol = 1L,
      dimnames = list(genes, "coef")
    ),
    p.value = matrix(
      pvalues,
      ncol = 1L,
      dimnames = list(genes, "pval")
    )
  )
}

make_collect_guideGIs <- function(data, space) {
  methods::new(
    "gRNA_GI",
    data = data,
    space = space,
    replicates = "replicate",
    block_layer = character(),
    blocks = character(),
    use_blocks = FALSE,
    block_description = character(),
    collapse = character()
  )
}

make_collect_screen <- function(
  class,
  guideGIs,
  limma_models,
  screen_attr = list(),
  metadata = list(),
  symmGeneGIs = data.table::data.table()
) {
  args <- list(
    Class = class,
    guideLFCs = methods::new(
      "gRNA_LFC",
      data = array(numeric(), dim = 0),
      space = character(),
      replicates = character()
    ),
    guideGIs = guideGIs,
    limma_models = limma_models,
    geneGIs = array(numeric(), dim = 0),
    screen_attr = screen_attr,
    dupCorrelation = numeric(),
    metadata = metadata,
    checks = list(),
    errors = list()
  )

  if (identical(class, "PosAgnMultiplexScreen")) {
    args$symmGeneGIs <- symmGeneGIs
  }

  do.call(methods::new, args)
}

make_fixed_pair_collect_data <- function() {
  matrix(
    seq_len(9L),
    nrow = 3L,
    dimnames = list(
      gene_pair = c("A;B", "A;C", "B;C"),
      replicate = c("r1", "r2", "r3")
    )
  )
}

make_multiplex_collect_data <- function() {
  array(
    seq_len(18L),
    dim = c(3L, 3L, 2L),
    dimnames = list(
      query_gene = c("A", "B", "C"),
      library_gene = c("A", "B", "C"),
      replicate = c("r1", "r2")
    )
  )
}

test_that("collect_GIs extracts fixed-pair coefficients, p-values, and adjusted FDR", {
  guide_data <- make_fixed_pair_collect_data()
  screen <- make_collect_screen(
    class = "FixedPairScreen",
    guideGIs = make_collect_guideGIs(guide_data, space = "gene_pair"),
    limma_models = make_collect_model(
      genes = rownames(guide_data),
      coefficients = c(0.5, -0.25, 1.2),
      pvalues = c(0.01, 0.20, 0.03)
    )
  )

  result <- collect_GIs(screen, FDR_method = "BH")

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(dim(result@geneGIs), c(3L, 3L))
  expect_equal(rownames(result@geneGIs), rownames(guide_data))
  expect_equal(colnames(result@geneGIs), c("GI", "pval", "FDR"))
  expect_equal(unname(result@geneGIs[, "GI"]), c(0.5, -0.25, 1.2))
  expect_equal(unname(result@geneGIs[, "pval"]), c(0.01, 0.20, 0.03))
  expect_equal(
    unname(result@geneGIs[, "FDR"]),
    stats::p.adjust(c(0.01, 0.20, 0.03), method = "BH")
  )
})

test_that("collect_GIs respects requested FDR method for fixed-pair screens", {
  guide_data <- make_fixed_pair_collect_data()
  screen <- make_collect_screen(
    class = "FixedPairScreen",
    guideGIs = make_collect_guideGIs(guide_data, space = "gene_pair"),
    limma_models = make_collect_model(
      genes = rownames(guide_data),
      coefficients = c(0.5, -0.25, 1.2),
      pvalues = c(0.01, 0.20, 0.03)
    )
  )

  result <- collect_GIs(screen, FDR_method = "bonferroni")

  expect_equal(
    unname(result@geneGIs[, "FDR"]),
    stats::p.adjust(c(0.01, 0.20, 0.03), method = "bonferroni")
  )
})

test_that("collect_GIs rejects unknown FDR methods", {
  guide_data <- make_fixed_pair_collect_data()
  screen <- make_collect_screen(
    class = "FixedPairScreen",
    guideGIs = make_collect_guideGIs(guide_data, space = "gene_pair"),
    limma_models = make_collect_model(
      rownames(guide_data),
      c(1, 2, 3),
      c(0.1, 0.2, 0.3)
    )
  )

  expect_error(
    collect_GIs(screen, FDR_method = "not-a-method"),
    "Unknown FDR method provided"
  )
})

test_that("collect_GIs builds a query-by-library-by-variable array for multiplex screens", {
  guide_data <- make_multiplex_collect_data()
  screen <- make_collect_screen(
    class = "MultiplexScreen",
    guideGIs = make_collect_guideGIs(
      guide_data,
      space = c("query_gene", "library_gene")
    ),
    limma_models = list(
      A = make_collect_model(
        c("A", "B", "C"),
        c(0.1, 0.2, 0.3),
        c(0.01, 0.02, 0.03)
      ),
      B = make_collect_model(
        c("A", "B", "C"),
        c(-0.1, -0.2, -0.3),
        c(0.04, 0.05, 0.06)
      ),
      C = make_collect_model(
        c("A", "B", "C"),
        c(0.7, 0.8, 0.9),
        c(0.07, 0.08, 0.09)
      )
    ),
    screen_attr = list(
      query_genes = c("A", "B", "C"),
      library_genes = c("A", "B", "C")
    )
  )

  result <- collect_GIs(screen, FDR_method = "BH")

  expect_s4_class(result, "MultiplexScreen")
  expect_equal(dim(result@geneGIs), c(3L, 3L, 3L))
  expect_equal(dimnames(result@geneGIs)[[1L]], c("A", "B", "C"))
  expect_equal(dimnames(result@geneGIs)[[2L]], c("A", "B", "C"))
  expect_equal(dimnames(result@geneGIs)[[3L]], c("GI", "pval", "FDR"))
  expect_equal(result@geneGIs["A", , "GI"], c(A = 0.1, B = 0.2, C = 0.3))
  expect_equal(result@geneGIs["B", , "GI"], c(A = -0.1, B = -0.2, C = -0.3))
  expect_equal(
    unname(result@geneGIs["A", , "FDR"]),
    stats::p.adjust(c(0.01, 0.02, 0.03), method = "BH")
  )
})

test_that("collect_GIs keeps failed multiplex models as NA rows", {
  guide_data <- make_multiplex_collect_data()
  screen <- make_collect_screen(
    class = "MultiplexScreen",
    guideGIs = make_collect_guideGIs(
      guide_data,
      space = c("query_gene", "library_gene")
    ),
    limma_models = list(
      A = make_collect_model(
        c("A", "B", "C"),
        c(0.1, 0.2, 0.3),
        c(0.01, 0.02, 0.03)
      ),
      B = NULL,
      C = make_collect_model(
        c("A", "B", "C"),
        c(0.7, 0.8, 0.9),
        c(0.07, 0.08, 0.09)
      )
    ),
    screen_attr = list(
      query_genes = c("A", "B", "C"),
      library_genes = c("A", "B", "C")
    )
  )

  result <- collect_GIs(screen)

  expect_true(all(is.na(result@geneGIs["B", , "GI"])))
  expect_true(all(is.na(result@geneGIs["B", , "pval"])))
  expect_true(all(is.na(result@geneGIs["B", , "FDR"])))
  expect_equal(result@geneGIs["A", "A", "GI"], 0.1)
})

test_that("collect_GIs warns when multiplex output loses guide-level genes", {
  guide_data <- make_multiplex_collect_data()
  screen <- make_collect_screen(
    class = "MultiplexScreen",
    guideGIs = make_collect_guideGIs(
      guide_data,
      space = c("query_gene", "library_gene")
    ),
    limma_models = list(
      A = make_collect_model(c("A", "B"), c(0.1, 0.2), c(0.01, 0.02)),
      B = make_collect_model(c("A", "B"), c(-0.1, -0.2), c(0.03, 0.04))
    ),
    screen_attr = list(
      query_genes = c("A", "B"),
      library_genes = c("A", "B")
    )
  )

  expect_warning(
    result <- collect_GIs(screen),
    "Some genes were lost"
  )
  expect_equal(dim(result@geneGIs), c(2L, 2L, 3L))
})

test_that("collect_GIs creates symmetrized output for position-agnostic multiplex screens", {
  guide_data <- make_multiplex_collect_data()
  screen <- make_collect_screen(
    class = "PosAgnMultiplexScreen",
    guideGIs = make_collect_guideGIs(
      guide_data,
      space = c("query_gene", "library_gene")
    ),
    limma_models = list(
      A = make_collect_model(
        c("A", "B", "C"),
        c(0.0, 0.2, 0.4),
        c(0.50, 0.01, 0.02)
      ),
      B = make_collect_model(
        c("A", "B", "C"),
        c(0.3, 0.0, 0.5),
        c(0.03, 0.50, 0.04)
      ),
      C = make_collect_model(
        c("A", "B", "C"),
        c(0.6, 0.7, 0.0),
        c(0.05, 0.06, 0.50)
      )
    ),
    screen_attr = list(
      query_genes = c("A", "B", "C"),
      library_genes = c("A", "B", "C"),
      unique_pairs = c("A;B", "A;C", "B;C")
    ),
    symmGeneGIs = data.table::data.table()
  )

  result <- collect_GIs(screen, FDR_method = "BH")

  expect_s4_class(result, "PosAgnMultiplexScreen")
  expect_s3_class(result@symmGeneGIs, "data.table")
  expect_named(
    result@symmGeneGIs,
    c("gene_pair", "query_gene", "library_gene", "GI", "GI_z", "pval", "FDR")
  )
  expect_equal(result@symmGeneGIs$gene_pair, c("A;B", "A;C", "B;C"))
  expect_equal(result@symmGeneGIs$query_gene, c("A", "A", "B"))
  expect_equal(result@symmGeneGIs$library_gene, c("B", "C", "C"))
  expect_equal(result@symmGeneGIs$GI, c(0.2, 0.4, 0.5))
  expect_equal(result@symmGeneGIs$pval, c(0.01, 0.02, 0.04))
  expect_equal(
    result@symmGeneGIs$FDR,
    balanced_FDR(
      pairs = c("A;B", "A;C", "B;C"),
      pval_array = result@geneGIs[,, "pval"],
      fdr_method = "BH"
    )
  )
  expect_equal(result@symmGeneGIs$GI_z, z_transform(c(0.2, 0.4, 0.5)))
})


test_that("collect_GIs applies one global FDR correction in canonical pair order", {
  pairs <- c("B;C", "A;B", "A;C")
  coefficients <- c(-0.5, 0.8, 0.1)
  pvalues <- c(0.001, 0.02, 0.40)
  guide_data <- matrix(
    seq_len(12L),
    nrow = 3L,
    dimnames = list(
      gene_pair = pairs,
      replicate = paste0("r", seq_len(4L))
    )
  )
  screen <- make_collect_screen(
    class = "PosAgnMultiplexScreen",
    guideGIs = make_collect_guideGIs(guide_data, space = "gene_pair"),
    limma_models = make_collect_model(
      genes = pairs,
      coefficients = coefficients,
      pvalues = pvalues
    ),
    screen_attr = list(unique_pairs = pairs),
    metadata = list(symmetric_analysis_method = "global_preaverage")
  )

  result <- collect_GIs(screen, FDR_method = "BH")
  expected_fdr <- stats::p.adjust(pvalues, method = "BH")

  expect_equal(result@symmGeneGIs$gene_pair, pairs)
  expect_equal(result@symmGeneGIs$query_gene, c("B", "A", "A"))
  expect_equal(result@symmGeneGIs$library_gene, c("C", "B", "C"))
  expect_equal(result@symmGeneGIs$GI, coefficients)
  expect_equal(result@symmGeneGIs$pval, pvalues)
  expect_equal(result@symmGeneGIs$FDR, expected_fdr)
  expect_equal(result@symmGeneGIs$GI_z, z_transform(coefficients))

  expect_equal(rownames(result@geneGIs), pairs)
  expect_equal(colnames(result@geneGIs), c("GI", "pval", "FDR"))
  expect_equal(unname(result@geneGIs[, "GI"]), coefficients)
  expect_equal(unname(result@geneGIs[, "pval"]), pvalues)
  expect_equal(unname(result@geneGIs[, "FDR"]), expected_fdr)
})


test_that("global position-agnostic collection respects the requested FDR method", {
  pairs <- c("A;B", "A;C", "B;C")
  pvalues <- c(0.01, 0.04, 0.20)
  guide_data <- matrix(
    seq_len(12L),
    nrow = 3L,
    dimnames = list(
      gene_pair = pairs,
      replicate = paste0("r", seq_len(4L))
    )
  )
  screen <- make_collect_screen(
    class = "PosAgnMultiplexScreen",
    guideGIs = make_collect_guideGIs(guide_data, space = "gene_pair"),
    limma_models = make_collect_model(
      genes = pairs,
      coefficients = c(0.1, 0.2, 0.3),
      pvalues = pvalues
    ),
    screen_attr = list(unique_pairs = pairs),
    metadata = list(symmetric_analysis_method = "global_preaverage")
  )

  result <- collect_GIs(screen, FDR_method = "bonferroni")

  expect_equal(
    result@symmGeneGIs$FDR,
    stats::p.adjust(pvalues, method = "bonferroni")
  )
})


test_that("global position-agnostic collection rejects reordered model pairs", {
  expected_pairs <- c("A;B", "A;C", "B;C")
  model_pairs <- c("A;C", "A;B", "B;C")
  guide_data <- matrix(
    seq_len(12L),
    nrow = 3L,
    dimnames = list(
      gene_pair = expected_pairs,
      replicate = paste0("r", seq_len(4L))
    )
  )
  screen <- make_collect_screen(
    class = "PosAgnMultiplexScreen",
    guideGIs = make_collect_guideGIs(guide_data, space = "gene_pair"),
    limma_models = make_collect_model(
      genes = model_pairs,
      coefficients = c(0.1, 0.2, 0.3),
      pvalues = c(0.01, 0.04, 0.20)
    ),
    screen_attr = list(unique_pairs = expected_pairs),
    metadata = list(symmetric_analysis_method = "global_preaverage")
  )

  expect_error(
    collect_GIs(screen),
    "global limma output rows do not match"
  )
})
