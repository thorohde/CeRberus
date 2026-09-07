make_gRNA_GI_for_fill <- function(
  space = "gene_pair",
  replicates = c("guide_pair", "bio_rep", "tech_rep"),
  block_layer = character(),
  blocks = character(),
  use_blocks = FALSE,
  block_description = character(),
  collapse = character()
) {
  methods::new(
    "gRNA_GI",
    data = array(numeric(), dim = 0),
    space = space,
    replicates = replicates,
    block_layer = block_layer,
    blocks = blocks,
    use_blocks = use_blocks,
    block_description = block_description,
    collapse = collapse
  )
}

make_fixed_pair_fill_scores <- function() {
  expand.grid(
    gene_pair = c("A;C", "B;D"),
    guide_pair = c("g1", "g2"),
    bio_rep = c("b1", "b2"),
    tech_rep = c("t1", "t2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    transform(
      GI = seq_len(16L),
      alternative_score = seq_len(16L) + 100L
    )
}

make_multiplex_fill_scores <- function() {
  expand.grid(
    query_gene = c("G1", "G2"),
    library_gene = c("G1", "G2", "G3"),
    guide_pair = c("g1", "g2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    transform(GI = seq_len(12L))
}

make_multiplex_fill_scores_with_all_na_gene_slices <- function() {
  input <- expand.grid(
    query_gene = c("G1", "G2", "G3"),
    library_gene = c("G1", "G2", "G3", "G4"),
    guide_pair = c("g1", "g2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  input$GI <- seq_len(nrow(input))
  input$GI[input$query_gene == "G3" | input$library_gene == "G4"] <- NA_real_

  input
}

test_that("fill_gRNA_GIs fills fixed-pair guide GI arrays from long input", {
  input <- make_fixed_pair_fill_scores()
  object <- make_gRNA_GI_for_fill()

  result <- CeRberus:::fill_gRNA_GIs(object, input)

  expect_equal(
    list(
      class = class(result)[[1]],
      dimensions = dim(result@data),
      dimnames = unname(dimnames(result@data)),
      first = result@data["A;C", "g1", "b1", "t1"],
      last = result@data["B;D", "g2", "b2", "t2"]
    ),
    list(
      class = "gRNA_GI",
      dimensions = c(2L, 2L, 2L, 2L),
      dimnames = list(
        c("A;C", "B;D"),
        c("g1", "g2"),
        c("b1", "b2"),
        c("t1", "t2")
      ),
      first = 1L,
      last = 16L
    )
  )
})

test_that("fill_gRNA_GIs fills multiplex arrays and retains all-NA genes", {
  input <- make_multiplex_fill_scores_with_all_na_gene_slices()
  object <- make_gRNA_GI_for_fill(
    space = c("query_gene", "library_gene"),
    replicates = "guide_pair"
  )

  result <- CeRberus:::fill_gRNA_GIs(object, input)

  expect_equal(
    list(
      dimensions = dim(result@data),
      dimnames = unname(dimnames(result@data)),
      first = result@data["G1", "G1", "g1"],
      query_all_na = all(is.na(result@data["G3", , ])),
      library_all_na = all(is.na(result@data[, "G4", ])),
      query_has_values = !all(is.na(result@data["G1", , ])),
      library_has_values = !all(is.na(result@data[, "G1", ]))
    ),
    list(
      dimensions = c(3L, 4L, 2L),
      dimnames = list(
        c("G1", "G2", "G3"),
        c("G1", "G2", "G3", "G4"),
        c("g1", "g2")
      ),
      first = 1,
      query_all_na = TRUE,
      library_all_na = TRUE,
      query_has_values = TRUE,
      library_has_values = TRUE
    )
  )
})

test_that("fill_gRNA_GIs supports custom values without changing metadata", {
  input <- make_fixed_pair_fill_scores()
  object <- make_gRNA_GI_for_fill(
    replicates = "guide_pair",
    block_layer = "guide_pair",
    blocks = "none",
    use_blocks = TRUE,
    block_description = "old_description",
    collapse = character()
  )
  input <- input[input$bio_rep == "b1" & input$tech_rep == "t1", ]

  result <- CeRberus:::fill_gRNA_GIs(
    object,
    input,
    value_var = "alternative_score"
  )

  expect_equal(
    list(
      values = unname(result@data[cbind(c(1L, 2L), c(1L, 2L))]),
      metadata = list(
        space = result@space,
        replicates = result@replicates,
        block_layer = result@block_layer,
        blocks = result@blocks,
        use_blocks = result@use_blocks,
        block_description = result@block_description,
        collapse = result@collapse
      )
    ),
    list(
      values = c(101, 104),
      metadata = list(
        space = object@space,
        replicates = object@replicates,
        block_layer = object@block_layer,
        blocks = object@blocks,
        use_blocks = object@use_blocks,
        block_description = object@block_description,
        collapse = object@collapse
      )
    )
  )
})

test_that("fill_gRNA_GIs represents missing input combinations as NA", {
  input <- make_fixed_pair_fill_scores()
  input <- input[
    !(input$gene_pair == "B;D" &
      input$guide_pair == "g2" &
      input$bio_rep == "b2" &
      input$tech_rep == "t2"),
  ]
  object <- make_gRNA_GI_for_fill()

  result <- CeRberus:::fill_gRNA_GIs(object, input)

  expect_equal(
    c(
      missing = result@data["B;D", "g2", "b2", "t2"],
      present = result@data["A;C", "g1", "b1", "t1"]
    ),
    c(missing = NA_real_, present = 1)
  )
})
