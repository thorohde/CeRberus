make_screen_for_attributes <- function(
  input,
  screen_attr = methods::new("ScreenDesign")
) {
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

test_that("get_screen_attributes replaces the design with derived attributes", {
  input <- make_attribute_input()
  screen <- make_screen_for_attributes(
    input,
    screen_attr = make_screen_design(query_genes = "existing")
  )

  result <- get_screen_attributes(screen)

  expect_equal(
    list(design = result@screen_attr, input = result@metadata$input),
    list(
      design = make_screen_design(
        query_genes = c("A", "B", "C", "D"),
        library_genes = c("B", "C", "A", "E"),
        all_pairs = c("A;B", "A;C", "B;A", "C;A", "C;E", "D;E"),
        observations_per_query = c(2L, 1L, 2L, 1L)
      ),
      input = input
    )
  )
})

test_that("get_screen_attributes handles duplicated rows without duplicating gene sets or pair lists", {
  input <- make_attribute_input()
  input <- rbind(input, input[1L])
  screen <- make_screen_for_attributes(input)

  result <- get_screen_attributes(screen)

  expect_equal(
    result@screen_attr,
    make_screen_design(
      query_genes = c("A", "B", "C", "D"),
      library_genes = c("B", "C", "A", "E"),
      all_pairs = c("A;B", "A;C", "B;A", "C;A", "C;E", "D;E"),
      observations_per_query = c(3L, 1L, 2L, 1L)
    )
  )
})
