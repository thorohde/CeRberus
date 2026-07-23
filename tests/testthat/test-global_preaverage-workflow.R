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

  result <- compute_dupCorrelation(result)
  result <- compute_models(result)
  result <- collect_GIs(result, FDR_method = "BH")
  output <- GI_df(result)

  expected_pairs <- result@screen_attr$unique_pairs
  model_pvalues <- as.numeric(result@limma_models$p.value[, 1L])

  expect_s4_class(result, "PosAgnMultiplexScreen")
  expect_identical(result@guideGIs@space, "gene_pair")
  expect_length(result@dupCorrelation, 1L)
  expect_true(inherits(result@limma_models, "MArrayLM"))
  expect_equal(rownames(result@limma_models$coefficients), expected_pairs)

  expect_s3_class(output, "data.table")
  expect_equal(output$gene_pair, expected_pairs)
  expect_equal(nrow(output), length(expected_pairs))
  expect_named(
    output,
    c("gene_pair", "query_gene", "library_gene", "GI", "GI_z", "pval", "FDR")
  )
  expect_equal(output$pval, model_pvalues)
  expect_equal(output$FDR, stats::p.adjust(model_pvalues, method = "BH"))
  expect_equal(output, result@symmGeneGIs)
})
