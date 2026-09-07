make_gRNA_GI_for_collapse <- function(data, collapse = character()) {
  methods::new(
    "gRNA_GI",
    data = data,
    space = c("query_gene", "library_gene"),
    replicates = c("bio_rep", "tech_rep"),
    block_layer = character(),
    blocks = character(),
    use_blocks = FALSE,
    block_description = character(),
    collapse = collapse
  )
}

make_collapse_array <- function() {
  array(
    1:16,
    dim = c(2, 2, 2, 2),
    dimnames = list(
      query_gene = c("A", "B"),
      library_gene = c("C", "D"),
      bio_rep = c("b1", "b2"),
      tech_rep = c("t1", "t2")
    )
  )
}

test_that("collapse_replicates returns the object unchanged when no collapse layer is set", {
  data <- make_collapse_array()
  object <- make_gRNA_GI_for_collapse(data)

  result <- collapse_replicates(object)

  expect_identical(result, object)
})

test_that("collapse_replicates averages selected replicate dimensions", {
  data <- make_collapse_array()
  cases <- list(
    tech_rep = list(
      collapse = "tech_rep",
      keep = c("query_gene", "library_gene", "bio_rep"),
      replicates = "bio_rep"
    ),
    bio_rep = list(
      collapse = "bio_rep",
      keep = c("query_gene", "library_gene", "tech_rep"),
      replicates = "tech_rep"
    ),
    both = list(
      collapse = c("bio_rep", "tech_rep"),
      keep = c("query_gene", "library_gene"),
      replicates = character()
    )
  )

  purrr::iwalk(cases, function(case, name) {
    result <- make_gRNA_GI_for_collapse(data, case$collapse) |>
      collapse_replicates()
    keep_margin <- match(case$keep, names(dimnames(data)))
    expected <- array(
      apply(data, keep_margin, mean, na.rm = TRUE),
      dim = vapply(dimnames(data)[case$keep], length, integer(1)),
      dimnames = dimnames(data)[case$keep]
    )

    expect_equal(
      list(data = result@data, replicates = result@replicates),
      list(data = expected, replicates = case$replicates),
      info = name
    )
  })
})

test_that("collapse_replicates handles missing values while averaging", {
  data <- make_collapse_array()
  data["A", "C", "b1", "t1"] <- NA_real_
  data["A", "C", "b1", "t2"] <- 10
  data["B", "D", "b2", ] <- NA_real_
  object <- make_gRNA_GI_for_collapse(data, collapse = "tech_rep")

  result <- collapse_replicates(object)

  expect_equal(
    c(
      partial = result@data["A", "C", "b1"],
      all_missing = result@data["B", "D", "b2"]
    ),
    c(partial = 10, all_missing = NaN)
  )
})

test_that("collapse_replicates rejects invalid collapse layers", {
  data <- make_collapse_array()
  cases <- list(
    non_replicate = list(
      action = function() {
        make_gRNA_GI_for_collapse(data, collapse = "query_gene")
      },
      error = "cannot be collapsed as replicates"
    ),
    unknown_layer = list(
      action = function() {
        make_gRNA_GI_for_collapse(
          data,
          collapse = c("tech_rep", "not_a_layer")
        ) |>
          collapse_replicates()
      },
      error = "not_a_layer"
    )
  )

  purrr::iwalk(cases, function(case, name) {
    expect_error(case$action(), case$error, info = name)
  })
})
