make_gRNA_GI_for_compute_models <- function(
  data,
  space,
  replicates = "replicate",
  blocks = c("b1", "b1", "b2", "b2")
) {
  methods::new(
    "gRNA_GI",
    data = data,
    space = space,
    replicates = replicates,
    block_layer = "replicate",
    blocks = blocks,
    use_blocks = TRUE,
    block_description = dimnames(data)[[length(space) + 1L]],
    collapse = character()
  )
}

make_compute_models_screen <- function(
  class,
  guideGIs,
  dupCorrelation,
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
    limma_models = list(),
    geneGIs = array(numeric(), dim = 0),
    screen_attr = screen_attr,
    dupCorrelation = dupCorrelation,
    metadata = metadata,
    checks = list(),
    errors = list()
  )

  if (identical(class, "PosAgnMultiplexScreen")) {
    args$symmGeneGIs <- symmGeneGIs
  }

  do.call(methods::new, args)
}

make_fixed_pair_model_matrix <- function() {
  matrix(
    c(
      1.0,
      1.2,
      1.8,
      2.1,
      2.0,
      2.2,
      2.9,
      3.2,
      3.0,
      3.1,
      3.8,
      4.0
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(
      gene_pair = c("A;C", "B;D", "E;F"),
      replicate = paste0("r", seq_len(4L))
    )
  )
}

make_multiplex_model_array <- function() {
  array(
    c(
      1.0,
      1.2,
      1.8,
      2.1,
      2.0,
      2.1,
      2.8,
      3.0,
      3.0,
      3.3,
      3.7,
      4.1,
      4.0,
      4.1,
      4.8,
      5.2,
      5.0,
      5.2,
      5.9,
      6.1,
      6.0,
      6.2,
      6.7,
      7.0
    ),
    dim = c(2L, 3L, 4L),
    dimnames = list(
      query_gene = c("Q1", "Q2"),
      library_gene = c("L1", "L2", "L3"),
      replicate = paste0("r", seq_len(4L))
    )
  )
}

test_that("compute_models fits and stores a limma model for fixed-pair screens", {
  data <- make_fixed_pair_model_matrix()
  guideGIs <- make_gRNA_GI_for_compute_models(
    data = data,
    space = "gene_pair"
  )
  screen <- make_compute_models_screen(
    class = "FixedPairScreen",
    guideGIs = guideGIs,
    dupCorrelation = 0.1
  )

  result <- compute_models(screen)

  expect_s4_class(result, "FixedPairScreen")
  expect_true(inherits(result@limma_models, "MArrayLM"))
  expect_equal(rownames(result@limma_models$coefficients), rownames(data))
  expect_equal(nrow(result@limma_models$coefficients), nrow(data))
  expect_equal(nrow(result@limma_models$p.value), nrow(data))
  expect_equal(result@guideGIs, guideGIs)
  expect_equal(result@dupCorrelation, 0.1)
})

test_that("compute_models fits one limma model per multiplex query gene", {
  data <- make_multiplex_model_array()
  guideGIs <- make_gRNA_GI_for_compute_models(
    data = data,
    space = c("query_gene", "library_gene")
  )
  screen <- make_compute_models_screen(
    class = "MultiplexScreen",
    guideGIs = guideGIs,
    dupCorrelation = c(Q1 = 0.05, Q2 = 0.1),
    screen_attr = list(
      query_genes = c("Q1", "Q2"),
      library_genes = c("L1", "L2", "L3")
    )
  )

  result <- compute_models(screen)

  expect_s4_class(result, "MultiplexScreen")
  expect_named(result@limma_models, c("Q1", "Q2"))
  expect_true(all(purrr::map_lgl(result@limma_models, inherits, "MArrayLM")))
  expect_equal(
    rownames(result@limma_models$Q1$coefficients),
    c("L1", "L2", "L3")
  )
  expect_equal(
    rownames(result@limma_models$Q2$coefficients),
    c("L1", "L2", "L3")
  )
  expect_equal(result@errors$query_genes_not_usable, character())
  expect_named(result@errors$GI_computation_errors, c("Q1", "Q2"))
  expect_true(all(purrr::map_lgl(result@errors$GI_computation_errors, is.null)))
})


test_that("compute_models fits one global model for global_preaverage screens", {
  data <- make_fixed_pair_model_matrix()
  guideGIs <- make_gRNA_GI_for_compute_models(
    data = data,
    space = "gene_pair"
  )
  screen <- make_compute_models_screen(
    class = "PosAgnMultiplexScreen",
    guideGIs = guideGIs,
    dupCorrelation = 0.1,
    screen_attr = list(unique_pairs = rownames(data)),
    metadata = list(symmetric_analysis_method = "global_preaverage")
  )

  result <- compute_models(screen)

  expect_s4_class(result, "PosAgnMultiplexScreen")
  expect_true(inherits(result@limma_models, "MArrayLM"))
  expect_equal(rownames(result@limma_models$coefficients), rownames(data))
  expect_equal(nrow(result@limma_models$coefficients), nrow(data))
  expect_equal(nrow(result@limma_models$p.value), nrow(data))
})


test_that("compute_models retains per-query models for preaverage screens", {
  data <- make_multiplex_model_array()
  guideGIs <- make_gRNA_GI_for_compute_models(
    data = data,
    space = c("query_gene", "library_gene")
  )
  screen <- make_compute_models_screen(
    class = "PosAgnMultiplexScreen",
    guideGIs = guideGIs,
    dupCorrelation = c(Q1 = 0.05, Q2 = 0.1),
    screen_attr = list(
      query_genes = c("Q1", "Q2"),
      library_genes = c("L1", "L2", "L3")
    ),
    metadata = list(symmetric_analysis_method = "preaverage")
  )

  result <- compute_models(screen)

  expect_s4_class(result, "PosAgnMultiplexScreen")
  expect_named(result@limma_models, c("Q1", "Q2"))
  expect_true(all(purrr::map_lgl(result@limma_models, inherits, "MArrayLM")))
})

test_that("compute_models stores failed multiplex query models and warns", {
  data <- make_multiplex_model_array()
  guideGIs <- make_gRNA_GI_for_compute_models(
    data = data,
    space = c("query_gene", "library_gene")
  )
  screen <- make_compute_models_screen(
    class = "MultiplexScreen",
    guideGIs = guideGIs,
    dupCorrelation = c(Q1 = 0.05, Q_missing = 0.1),
    screen_attr = list(
      query_genes = c("Q1", "Q_missing"),
      library_genes = c("L1", "L2", "L3")
    )
  )

  expect_warning(
    result <- compute_models(screen),
    "Failed computing GIs for 1 genes"
  )

  expect_true(inherits(result@limma_models$Q1, "MArrayLM"))
  expect_null(result@limma_models$Q_missing)
  expect_equal(result@errors$query_genes_not_usable, "Q_missing")
  expect_null(result@errors$GI_computation_errors$Q1)
  expect_s3_class(result@errors$GI_computation_errors$Q_missing, "error")
})

test_that("compute_models errors directly for invalid fixed-pair model input", {
  data <- make_fixed_pair_model_matrix()
  guideGIs <- make_gRNA_GI_for_compute_models(
    data = data,
    space = "gene_pair",
    blocks = c("b1", "b2")
  )
  screen <- make_compute_models_screen(
    class = "FixedPairScreen",
    guideGIs = guideGIs,
    dupCorrelation = 0.1
  )

  expect_error(compute_models(screen))
})
