test_that("sort_gene_pairs supports vector and concatenated-pair inputs", {
  g1 <- c("RB1", "NOTCH1", "TTN", "MSH2", "FANCD2")
  g2 <- c("RNF43", "NEK1", "HLTF", "REV3L", "PAPD7")
  cases <- list(
    separate_vectors = list(
      action = function() sort_gene_pairs(g1, g2),
      expected = c(
        "RB1;RNF43",
        "NEK1;NOTCH1",
        "HLTF;TTN",
        "MSH2;REV3L",
        "FANCD2;PAPD7"
      )
    ),
    inverted = list(
      action = function() sort_gene_pairs(g1, g2, invert = TRUE),
      expected = c(
        "RNF43;RB1",
        "NOTCH1;NEK1",
        "TTN;HLTF",
        "REV3L;MSH2",
        "PAPD7;FANCD2"
      )
    ),
    output_separator = list(
      action = function() sort_gene_pairs(c("B", "A"), c("A", "C"), sep = "|"),
      expected = c("A|B", "A|C")
    ),
    concatenated_pairs = list(
      action = function() {
        sort_gene_pairs(pairs = c("RB1;RNF43", "NOTCH1;NEK1", "TTN;HLTF"))
      },
      expected = c("RB1;RNF43", "NEK1;NOTCH1", "HLTF;TTN")
    ),
    input_separator = list(
      action = function() {
        sort_gene_pairs(pairs = c("B|A", "A|C"), pair_sep = "\\|", sep = ";")
      },
      expected = c("A;B", "A;C")
    ),
    explicit_vectors = list(
      action = function() {
        sort_gene_pairs(
          g1 = c("B", "D"),
          g2 = c("A", "C"),
          pairs = c("X;Y", "Y;Z")
        )
      },
      expected = c("A;B", "C;D")
    ),
    duplicate_and_self_pairs = list(
      action = function() {
        sort_gene_pairs(c("B", "A", "A"), c("A", "B", "A"))
      },
      expected = c("A;B", "A;B", "A;A")
    )
  )

  purrr::iwalk(cases, function(case, name) {
    expect_equal(case$action(), case$expected, info = name)
  })
})

test_that("sort_gene_pairs validates required inputs and options", {
  cases <- list(
    absent = list(
      action = function() sort_gene_pairs(),
      error = "No input pairs given"
    ),
    missing_g2 = list(
      action = function() sort_gene_pairs(g1 = "A"),
      error = "Either two gene vectors or a gene pair vector required"
    ),
    missing_g1 = list(
      action = function() sort_gene_pairs(g2 = "B"),
      error = "Either two gene vectors or a gene pair vector required"
    ),
    unequal_lengths = list(
      action = function() sort_gene_pairs(g1 = c("A", "B"), g2 = "C"),
      error = "g1 and g2 have to be of equal length"
    ),
    invert = list(
      action = function() sort_gene_pairs("A", "B", invert = "FALSE"),
      error = "invert needs to be logical"
    )
  )

  purrr::iwalk(cases, function(case, name) {
    expect_error(case$action(), case$error, info = name)
  })
})
