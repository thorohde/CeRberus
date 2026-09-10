make_global_preaverage_workflow_scores <- function(n_genes = 20L) {
  genes <- paste0("G", seq_len(n_genes))
  input <- expand.grid(
    query_gene = genes,
    library_gene = genes,
    guide_pair = c("g1", "g2"),
    bio_rep = c("b1", "b2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  query_index <- match(input$query_gene, genes)
  library_index <- match(input$library_gene, genes)
  guide_offset <- ifelse(input$guide_pair == "g1", -0.15, 0.15)
  bio_offset <- ifelse(input$bio_rep == "b1", -0.05, 0.05)

  input$GI <- sin(query_index / 3) +
    cos(library_index / 4) +
    query_index * library_index / 100 +
    guide_offset +
    bio_offset

  input
}


test_that("global_preaverage runs end to end with one ordered pair output", {
  result <- GIScores(
    make_global_preaverage_workflow_scores(),
    pos_agnostic = TRUE,
    symmetric_analysis_method = "global_preaverage",
    block_layer = "guide_pair"
  )

  result <- compute_dup_correlation(result)
  result <- compute_models(result)
  result <- collect_gis(result, fdr_method = "BH")
  output <- gi_df(result)

  expected_pairs <- result@screen_attr$unique_pairs
  model_pvalues <- as.numeric(result@aggregatedLimmaModels$p.value[, 1L])

  expect_s4_class(result, "PosAgnMultiplexScreen")
  expect_identical(result@guideGIs@space, c("query_gene", "library_gene"))
  expect_identical(result@aggregatedGuideGIs@space, "gene_pair")
  expect_length(result@dupCorrelation, length(result@screen_attr$query_genes))
  expect_length(result@metadata$aggregated_dupCorrelation, 1L)
  expect_length(result@limma_models, 0L)
  expect_true(inherits(result@aggregatedLimmaModels, "MArrayLM"))
  expect_equal(
    rownames(result@aggregatedLimmaModels$coefficients),
    expected_pairs
  )

  expect_s3_class(output, "data.table")
  expect_equal(output$gene_pair, expected_pairs)
  expect_equal(nrow(output), length(expected_pairs))
  expect_named(
    output,
    c(
      "gene_pair",
      "query_gene",
      "library_gene",
      "GI",
      "GI_z",
      "pval",
      "FDR",
      "has_NTC"
    )
  )
  expect_equal(output$pval, model_pvalues)
  expect_equal(output$FDR, stats::p.adjust(model_pvalues, method = "BH"))
  expect_equal(output, append_ntc_pair_annotations(result@symmGeneGIs, result))
})

test_that("global_preaverage retains directional and aggregate results on request", {
  result <- GIScores(
    make_global_preaverage_workflow_scores(),
    pos_agnostic = TRUE,
    symmetric_analysis_method = "global_preaverage",
    retain_directional = TRUE,
    block_layer = "guide_pair"
  )

  result <- compute_dup_correlation(result)
  result <- compute_models(result)
  result <- collect_gis(result, fdr_method = "BH")
  output <- gi_df(result)

  expect_named(
    result@limma_models,
    result@screen_attr$query_genes
  )
  expect_true(all(purrr::map_lgl(result@limma_models, inherits, "MArrayLM")))
  expect_true(inherits(result@aggregatedLimmaModels, "MArrayLM"))
  expect_equal(
    dim(result@geneGIs),
    c(
      length(result@screen_attr$query_genes),
      length(result@screen_attr$library_genes),
      3L
    )
  )
  expect_true(all(c(
    "GI_ab",
    "pval_ab",
    "FDR_ab",
    "GI_ba",
    "pval_ba",
    "FDR_ba",
    "GI_aggregated",
    "pval_aggregated",
    "FDR_aggregated"
  ) %in% names(output)))
  expect_equal(output$GI, output$GI_aggregated)
  expect_equal(output$pval, output$pval_aggregated)
  expect_equal(output$FDR, output$FDR_aggregated)
  expect_true(result@metadata$multiple_testing$directional$retained)
  expect_identical(
    result@metadata$multiple_testing$aggregated$scope,
    "global_unordered_pairs"
  )
})
