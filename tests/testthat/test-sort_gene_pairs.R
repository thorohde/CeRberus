test_that("sort_gene_pairs sorts separate gene vectors alphabetically", {
  g1 <- c("RB1", "NOTCH1", "TTN", "MSH2", "FANCD2")
  g2 <- c("RNF43", "NEK1", "HLTF", "REV3L", "PAPD7")

  result <- sort_gene_pairs(g1, g2)

  expect_equal(
    result,
    c("RB1;RNF43", "NEK1;NOTCH1", "HLTF;TTN", "MSH2;REV3L", "FANCD2;PAPD7")
  )
})

test_that("sort_gene_pairs can invert alphabetical sorting", {
  g1 <- c("RB1", "NOTCH1", "TTN", "MSH2", "FANCD2")
  g2 <- c("RNF43", "NEK1", "HLTF", "REV3L", "PAPD7")

  result <- sort_gene_pairs(g1, g2, invert = TRUE)

  expect_equal(
    result,
    c("RNF43;RB1", "NOTCH1;NEK1", "TTN;HLTF", "REV3L;MSH2", "PAPD7;FANCD2")
  )
})

test_that("sort_gene_pairs supports a custom output separator", {
  result <- sort_gene_pairs(
    g1 = c("B", "A"),
    g2 = c("A", "C"),
    sep = "|"
  )

  expect_equal(result, c("A|B", "A|C"))
})

test_that("sort_gene_pairs sorts pre-concatenated pair vectors", {
  pairs <- c("RB1;RNF43", "NOTCH1;NEK1", "TTN;HLTF")

  result <- sort_gene_pairs(pairs = pairs)

  expect_equal(result, c("RB1;RNF43", "NEK1;NOTCH1", "HLTF;TTN"))
})

test_that("sort_gene_pairs supports a custom input pair separator", {
  pairs <- c("B|A", "A|C")

  result <- sort_gene_pairs(pairs = pairs, pair_sep = "\\|", sep = ";")

  expect_equal(result, c("A;B", "A;C"))
})

test_that("sort_gene_pairs prefers explicit gene vectors when both input styles are provided", {
  result <- sort_gene_pairs(
    g1 = c("B", "D"),
    g2 = c("A", "C"),
    pairs = c("X;Y", "Y;Z")
  )

  expect_equal(result, c("A;B", "C;D"))
})

test_that("sort_gene_pairs handles duplicated and identical gene names", {
  result <- sort_gene_pairs(
    g1 = c("B", "A", "A"),
    g2 = c("A", "B", "A")
  )

  expect_equal(result, c("A;B", "A;B", "A;A"))
})

test_that("sort_gene_pairs validates required input", {
  expect_error(
    sort_gene_pairs(),
    "No input pairs given"
  )

  expect_error(
    sort_gene_pairs(g1 = "A"),
    "Either two gene vectors or a gene pair vector required"
  )

  expect_error(
    sort_gene_pairs(g2 = "B"),
    "Either two gene vectors or a gene pair vector required"
  )
})

test_that("sort_gene_pairs validates vector lengths", {
  expect_error(
    sort_gene_pairs(g1 = c("A", "B"), g2 = "C"),
    "g1 and g2 have to be of equal length"
  )
})

test_that("sort_gene_pairs validates invert", {
  expect_error(
    sort_gene_pairs(g1 = "A", g2 = "B", invert = "FALSE"),
    "invert needs to be logical"
  )
})
