make_screen_for_attributes <- function(input, screen_attr = list(existing = "remove me")) {
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
    screen_attr = screen_attr,
    dupCorrelation = numeric(),
    metadata = list(input = input),
    checks = list(),
    errors = list()
  )
}

make_attribute_input <- function() {
  data.table::data.table(
    query_gene = c("A", "A", "B", "C", "C", "D"),
    library_gene = c("B", "C", "A", "A", "E", "E"),
    gene_pair = c("A;B", "A;C", "B;A", "C;A", "C;E", "D;E"),
    GI = seq_len(6L)
  )
}

test_that("get_screen_attributes computes gene sets and counts", {
  input <- make_attribute_input()
  screen <- make_screen_for_attributes(input)

  result <- get_screen_attributes(screen)

  expect_s4_class(result, "ScreenBase")
  expect_equal(result@screen_attr$contrasts, NULL)
  expect_equal(result@screen_attr$query_genes, c("A", "B", "C", "D"))
  expect_equal(result@screen_attr$library_genes, c("B", "C", "A", "E"))
  expect_equal(result@screen_attr$all_genes, c("A", "B", "C", "D", "E"))
  expect_equal(result@screen_attr$n_query_genes, 4L)
  expect_equal(result@screen_attr$n_lib_genes, 4L)
  expect_equal(result@screen_attr$n_all_genes, 5L)
})

test_that("get_screen_attributes reports genes missing from the opposite axis", {
  result <- get_screen_attributes(make_screen_for_attributes(make_attribute_input()))

  expect_equal(result@screen_attr$query_genes_not_in_lib, "D")
  expect_equal(result@screen_attr$library_genes_not_in_query, "E")
})

test_that("get_screen_attributes counts observations per query gene", {
  result <- get_screen_attributes(make_screen_for_attributes(make_attribute_input()))

  expect_equal(
    result@screen_attr$observations_per_query,
    c(A = 2L, B = 1L, C = 2L, D = 1L)
  )
})

test_that("get_screen_attributes stores directional and unordered gene pairs", {
  result <- get_screen_attributes(make_screen_for_attributes(make_attribute_input()))

  expect_equal(
    result@screen_attr$all_pairs,
    c("A;B", "A;C", "B;A", "C;A", "C;E", "D;E")
  )
  expect_equal(
    result@screen_attr$unique_pairs,
    c("A;B", "A;C", "C;E", "D;E")
  )
})

test_that("get_screen_attributes preserves metadata input and replaces previous attributes", {
  input <- make_attribute_input()
  screen <- make_screen_for_attributes(input)

  result <- get_screen_attributes(screen)

  expect_equal(result@metadata$input, input)
  expect_false("existing" %in% names(result@screen_attr))
})

test_that("get_screen_attributes handles duplicated rows without duplicating gene sets or pair lists", {
  input <- make_attribute_input()
  input <- rbind(input, input[1L])
  screen <- make_screen_for_attributes(input)

  result <- get_screen_attributes(screen)

  expect_equal(result@screen_attr$query_genes, c("A", "B", "C", "D"))
  expect_equal(result@screen_attr$library_genes, c("B", "C", "A", "E"))
  expect_equal(result@screen_attr$all_pairs, c("A;B", "A;C", "B;A", "C;A", "C;E", "D;E"))
  expect_equal(result@screen_attr$observations_per_query, c(A = 3L, B = 1L, C = 2L, D = 1L))
})

test_that("get_screen_attributes requires data.table input columns used by the method", {
  input <- as.data.frame(make_attribute_input())
  screen <- make_screen_for_attributes(input)

  expect_error(get_screen_attributes(screen))

  missing_gene_pair <- data.table::copy(make_attribute_input())
  missing_gene_pair[, gene_pair := NULL]
  screen <- make_screen_for_attributes(missing_gene_pair)

  expect_error(get_screen_attributes(screen))
})
