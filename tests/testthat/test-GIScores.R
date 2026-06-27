make_fixed_pair_scores <- function() {
  data.frame(
    query_gene = rep(c("A", "B"), each = 4),
    library_gene = rep(c("C", "D"), each = 4),
    bio_rep = rep(c("b1", "b2"), 4),
    tech_rep = rep(c("t1", "t2"), each = 2, times = 2),
    guide_pair = rep(c("g1", "g2"), 4),
    GI = seq_len(8)
  )
}

make_multiplex_scores <- function(n_genes = 20L) {
  genes <- paste0("G", seq_len(n_genes))
  input <- expand.grid(
    query_gene = genes,
    library_gene = genes,
    guide_pair = "g1",
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  input$GI <- seq_len(nrow(input))
  input
}

test_that("GIScores constructs a fixed-pair screen from guide-level scores", {
  input <- make_fixed_pair_scores()
  result <- GIScores(input, block_layer = "guide_pair")

  expect_s4_class(result, "FixedPairScreen")
  expect_s4_class(result@guideGIs, "gRNA_GI")
  expect_equal(result@guideGIs@space, "gene_pair")
  expect_equal(result@guideGIs@replicates, c("guide_pair", "tech_rep", "bio_rep"))
  expect_setequal(rownames(result@guideGIs@data), c("A;C", "B;D"))
  expect_equal(result@screen_attr$n_query_genes, 2L)
  expect_equal(result@screen_attr$n_lib_genes, 2L)
})

test_that("GIScores accepts data.tables without modifying the input", {
  input <- data.table::as.data.table(make_fixed_pair_scores())
  original <- data.table::copy(input)
  result <- GIScores(input, block_layer = "guide_pair")

  expect_equal(input, original)
  expect_true(data.table::is.data.table(result@metadata$input))
  expect_true("gene_pair" %in% names(result@metadata$input))
})

test_that("GIScores standardizes custom input column names", {
  input <- make_fixed_pair_scores()
  names(input) <- c("query", "library", "biological", "technical", "guide", "score")

  result <- GIScores(
    input,
    query_col = "query",
    lib_col = "library",
    bio_rep_col = "biological",
    tech_rep_col = "technical",
    guide_col = "guide",
    gi_col = "score",
    block_layer = "guide_pair"
  )

  expect_true(all(
    c("query_gene", "library_gene", "bio_rep", "tech_rep", "guide_pair", "GI") %in%
      names(result@metadata$input)
  ))
  expect_equal(result@metadata$input$GI, seq_len(8))
})

test_that("GIScores configures blocking from a retained replicate layer", {
  result <- GIScores(make_fixed_pair_scores(), block_layer = "guide_pair")

  expect_true(result@guideGIs@use_blocks)
  expect_length(result@guideGIs@blocks, ncol(result@guideGIs@data))
  expect_setequal(unique(result@guideGIs@blocks), c("g1", "g2"))
})

test_that("GIScores infers a multiplex screen for a complete gene grid", {
  result <- GIScores(make_multiplex_scores(), block_layer = "guide_pair")

  expect_s4_class(result, "MultiplexScreen")
  expect_equal(result@guideGIs@space, c("query_gene", "library_gene"))
  expect_equal(dim(result@guideGIs@data)[1:2], c(20L, 20L))
  expect_equal(result@screen_attr$n_all_genes, 20L)
  expect_true(result@checks$gene_sets_equal)
  expect_true(result@checks$sufficient_tests_per_query)
})

test_that("GIScores can force a multiplex-shaped input to fixed-pair mode", {
  expect_warning(
    result <- GIScores(
      make_multiplex_scores(),
      force_fixed_pair = TRUE,
      block_layer = "guide_pair"
    ),
    "Set up to use fixed pair structure"
  )

  expect_s4_class(result, "FixedPairScreen")
  expect_equal(result@guideGIs@space, "gene_pair")
})

test_that("GIScores creates a position-agnostic symmetric multiplex screen", {
  input <- make_multiplex_scores()
  input$GI <- as.numeric(factor(input$query_gene)) * 100 +
    as.numeric(factor(input$library_gene))
  result <- GIScores(
    input,
    pos_agnostic = TRUE,
    block_layer = "guide_pair"
  )

  expect_s4_class(result, "PosAgnMultiplexScreen")
  for (replicate_name in result@guideGIs@block_description) {
    expect_true(isSymmetric(result@guideGIs@data[, , replicate_name]))
  }
  expected <- mean(c(
    input$GI[input$query_gene == "G1" & input$library_gene == "G2"],
    input$GI[input$query_gene == "G2" & input$library_gene == "G1"]
  ))
  expect_equal(result@guideGIs@data["G1", "G2", 1], expected)
})

test_that("GIScores validates its input and required gene columns", {
  expect_error(
    GIScores(list(query_gene = "A", library_gene = "B")),
    "needs to be a data frame"
  )
  expect_error(
    GIScores(data.frame(library_gene = "B", GI = 1)),
    "query gene column"
  )
  expect_error(
    GIScores(data.frame(query_gene = "A", GI = 1)),
    "library gene column"
  )
})

test_that("GIScores rejects collapse layers that are not replicate dimensions", {
  expect_error(
    GIScores(make_fixed_pair_scores(), collapse_layers = "not_a_layer"),
    "not_a_layer"
  )
})
