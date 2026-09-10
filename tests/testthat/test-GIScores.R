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

test_that("GIScores constructs a fixed-pair screen with its inferred design", {
  input <- make_fixed_pair_scores()
  result <- GIScores(input, block_layer = "guide_pair")

  expect_s4_class(result, "FixedPairScreen")
  expect_s4_class(result@guideGIs, "gRNA_GI")
  expect_s4_class(result@guideLFCs, "gRNA_LFC")
  expect_s4_class(result@screen_attr, "ScreenDesign")
  expect_equal(result@guideGIs@space, "gene_pair")
  expect_equal(
    result@guideGIs@replicates,
    c("guide_pair", "tech_rep", "bio_rep")
  )
  expect_true(result@guideGIs@use_blocks)
  expect_length(result@guideGIs@blocks, ncol(result@guideGIs@data))
  expect_setequal(unique(result@guideGIs@blocks), c("g1", "g2"))
  expect_setequal(rownames(result@guideGIs@data), c("A;C", "B;D"))
  expect_length(result@guideLFCs@data, 0L)
  expect_equal(result@screen_attr$n_query_genes, 2L)
  expect_equal(result@screen_attr$n_lib_genes, 2L)
  expect_equal(result@screen_attr$query_genes, c("A", "B"))
  expect_equal(result@screen_attr$library_genes, c("C", "D"))
  expect_equal(result@screen_attr[["all_genes"]], c("A", "B", "C", "D"))
  expect_equal(
    result@screen_attr$observations_per_query,
    c(A = 4L, B = 4L)
  )
  expect_setequal(result@screen_attr$unique_pairs, c("A;C", "B;D"))
  expect_error(
    screen_attr(result) <- list(query_genes = "A"),
    "ScreenDesign"
  )
})

test_that("GIScores keeps optional guide LFCs aligned with guide GIs", {
  input <- make_fixed_pair_scores()
  input$LFC <- input$GI / 10

  results <- list(
    default = GIScores(input, block_layer = "guide_pair"),
    collapsed = GIScores(
      input,
      collapse_layers = "tech_rep",
      block_layer = "guide_pair"
    )
  )

  purrr::iwalk(results, function(result, strategy) {
    expect_identical(
      dim(result@guideLFCs@data),
      dim(result@guideGIs@data),
      info = strategy
    )
    expect_identical(
      dimnames(result@guideLFCs@data),
      dimnames(result@guideGIs@data),
      info = strategy
    )
    expect_identical(
      result@guideLFCs@space,
      result@guideGIs@space,
      info = strategy
    )
    expect_identical(
      result@guideLFCs@replicates,
      result@guideGIs@replicates,
      info = strategy
    )
  })

  result <- results$default
  expect_equal(
    result@guideLFCs@data["A;C", "g1_t1_b1"],
    input$LFC[
      input$query_gene == "A" &
        input$guide_pair == "g1" &
        input$bio_rep == "b1" &
        input$tech_rep == "t1"
    ]
  )
  expect_equal(
    result@guideLFCs@query_main_effects,
    c(
      A = mean(input$LFC[input$query_gene == "A"]),
      B = mean(input$LFC[input$query_gene == "B"])
    )
  )
  expect_equal(
    result@guideLFCs@library_main_effects,
    c(
      C = mean(input$LFC[input$library_gene == "C"]),
      D = mean(input$LFC[input$library_gene == "D"])
    )
  )
})

test_that("GIScores automatically computes multiplex gene main effects", {
  input <- make_multiplex_scores()
  input$LFC <- as.numeric(factor(input$query_gene)) +
    10 * as.numeric(factor(input$library_gene))

  result <- GIScores(
    input,
    screen_type = "multiplex",
    block_layer = "guide_pair"
  )

  query_genes <- unique(input$query_gene)
  expected_query <- stats::setNames(
    vapply(
      query_genes,
      \(.gene) mean(input$LFC[input$query_gene == .gene]),
      numeric(1L)
    ),
    query_genes
  )
  library_genes <- unique(input$library_gene)
  expected_library <- stats::setNames(
    vapply(
      library_genes,
      \(.gene) mean(input$LFC[input$library_gene == .gene]),
      numeric(1L)
    ),
    library_genes
  )

  expect_equal(result@guideLFCs@query_main_effects, expected_query)
  expect_equal(result@guideLFCs@library_main_effects, expected_library)
})

test_that("GIScores standardizes custom input columns including guide LFC", {
  input <- make_fixed_pair_scores()
  input$LFC <- input$GI / 10
  names(input) <- c(
    "query",
    "library",
    "biological",
    "technical",
    "guide",
    "score",
    "guide_log_fold_change"
  )
  input <- data.table::as.data.table(input)
  original <- data.table::copy(input)

  result <- GIScores(
    input,
    query_col = "query",
    lib_col = "library",
    bio_rep_col = "biological",
    tech_rep_col = "technical",
    guide_col = "guide",
    gi_col = "score",
    block_layer = "guide_pair",
    lfc_col = "guide_log_fold_change"
  )

  expect_equal(input, original)
  expect_true(data.table::is.data.table(result@metadata$input))
  expect_true(all(
    c(
      "query_gene",
      "library_gene",
      "bio_rep",
      "tech_rep",
      "guide_pair",
      "GI",
      "LFC",
      "gene_pair"
    ) %in%
      names(result@metadata$input)
  ))
  expect_equal(result@metadata$input$GI, seq_len(8))
  expect_equal(result@metadata$input$LFC, input$guide_log_fold_change)
  expect_identical(dim(result@guideLFCs@data), dim(result@guideGIs@data))
})

test_that("GIScores validates constructor input columns", {
  cases <- list(
    input_type = list(
      call = function() {
        GIScores(list(query_gene = "A", library_gene = "B"))
      },
      error = "needs to be a data frame"
    ),
    query_gene = list(
      call = function() GIScores(data.frame(library_gene = "B", GI = 1)),
      error = "query gene column"
    ),
    library_gene = list(
      call = function() GIScores(data.frame(query_gene = "A", GI = 1)),
      error = "library gene column"
    ),
    missing_lfc = list(
      call = function() {
        GIScores(make_fixed_pair_scores(), lfc_col = "missing_lfc")
      },
      error = "guide LFC column"
    ),
    non_numeric_lfc = list(
      call = function() {
        input <- make_fixed_pair_scores()
        input$LFC <- as.character(input$GI)
        GIScores(input)
      },
      error = "guide LFC column must be numeric"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_error(case$call(), case$error, info = case_name)
  })
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

test_that("GIScores supports explicit screen-type overrides", {
  cases <- list(
    fixed_pair = list(
      input = make_multiplex_scores(),
      inferred = "multiplex",
      class = "FixedPairScreen",
      space = "gene_pair"
    ),
    multiplex = list(
      input = make_fixed_pair_scores(),
      inferred = "fixed_pair",
      class = "MultiplexScreen",
      space = c("query_gene", "library_gene")
    )
  )

  purrr::iwalk(cases, function(case, screen_type) {
    expect_warning(
      result <- GIScores(
        case$input,
        screen_type = screen_type,
        block_layer = "guide_pair"
      ),
      paste0("overrides inferred screen type '", case$inferred, "'"),
      info = screen_type
    )

    expect_s4_class(result, case$class)
    expect_equal(result@guideGIs@space, case$space, info = screen_type)
    expect_identical(
      result@metadata$requested_screen_type,
      screen_type,
      info = screen_type
    )
    expect_identical(
      result@metadata$inferred_screen_type,
      case$inferred,
      info = screen_type
    )
    expect_identical(
      result@metadata$selected_screen_type,
      screen_type,
      info = screen_type
    )
  })
})

test_that("GIScores supports the deprecated force_fixed_pair override", {
  warning_messages <- character()
  result <- withCallingHandlers(
    GIScores(
      make_multiplex_scores(),
      force_fixed_pair = TRUE,
      block_layer = "guide_pair"
    ),
    warning = function(warning_condition) {
      warning_messages <<- c(
        warning_messages,
        conditionMessage(warning_condition)
      )
      invokeRestart("muffleWarning")
    }
  )

  expect_s4_class(result, "FixedPairScreen")
  expect_identical(result@metadata$requested_screen_type, "fixed_pair")
  expect_true(any(grepl(
    "force_fixed_pair is deprecated",
    warning_messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "overrides inferred screen type 'multiplex'",
    warning_messages,
    fixed = TRUE
  )))
})

test_that("GIScores validates constructor options", {
  cases <- list(
    screen_type = list(
      call = function() {
        GIScores(make_fixed_pair_scores(), screen_type = "unknown")
      },
      error = "should be one of"
    ),
    conflicting_screen_type = list(
      call = function() {
        suppressWarnings(GIScores(
          make_fixed_pair_scores(),
          force_fixed_pair = TRUE,
          screen_type = "multiplex"
        ))
      },
      error = "Do not use force_fixed_pair"
    ),
    collapse_layer = list(
      call = function() {
        GIScores(make_fixed_pair_scores(), collapse_layers = "not_a_layer")
      },
      error = "not_a_layer"
    ),
    retain_directional = list(
      call = function() {
        GIScores(make_fixed_pair_scores(), retain_directional = "yes")
      },
      error = "retain_directional must be TRUE or FALSE"
    )
  )

  purrr::iwalk(cases, function(case, case_name) {
    expect_error(case$call(), case$error, info = case_name)
  })
})

test_that("GIScores preserves the legacy positional argument order", {
  result <- GIScores(
    make_fixed_pair_scores(),
    "query_gene",
    "library_gene",
    "bio_rep",
    "tech_rep",
    "guide_pair",
    "GI",
    NULL,
    "guide_pair",
    FALSE,
    FALSE,
    "preaverage",
    FALSE
  )

  expect_s4_class(result, "FixedPairScreen")
  expect_identical(result@metadata$requested_screen_type, "auto")
  expect_false(result@metadata$retain_directional)
})

test_that("GIScores stores non-targeting controls in metadata", {
  result <- GIScores(
    make_fixed_pair_scores(),
    non_targeting_controls = c("NTC_1", "NTC_2")
  )

  expect_identical(
    result@metadata$non_targeting_controls,
    c("NTC_1", "NTC_2")
  )
})

test_that("GIScores validates non-targeting controls", {
  for (controls in list(character(), c("NTC", "NTC"), NA_character_, 1)) {
    expect_error(
      GIScores(
        make_fixed_pair_scores(),
        non_targeting_controls = controls
      ),
      "non_targeting_controls must be NULL or a non-empty character vector"
    )
  }
})

test_that("GIScores records unmatched non-targeting controls", {
  input <- make_fixed_pair_scores()
  input$query_gene[input$query_gene == "A"] <- "NTC"

  matched <- GIScores(input, non_targeting_controls = "NTC")
  expect_identical(matched@metadata$unmatched_non_targeting_controls, character())

  expect_warning(
    unmatched <- GIScores(
      input,
      non_targeting_controls = c("NTC", "NOT_IN_SCREEN")
    ),
    "Configured non-targeting controls not found"
  )
  expect_identical(
    unmatched@metadata$unmatched_non_targeting_controls,
    "NOT_IN_SCREEN"
  )
})

test_that("GIScores creates a position-agnostic symmetric multiplex screen", {
  input <- make_multiplex_scores()
  input$GI <- as.numeric(factor(input$query_gene)) *
    100 +
    as.numeric(factor(input$library_gene))
  input$LFC <- input$GI / 10
  result <- GIScores(
    input,
    pos_agnostic = TRUE,
    block_layer = "guide_pair"
  )

  expect_s4_class(result, "PosAgnMultiplexScreen")
  expect_false(result@metadata$retain_directional)
  expect_s4_class(result@aggregatedGuideGIs, "gRNA_GI")
  for (replicate_name in result@aggregatedGuideGIs@block_description) {
    expect_true(isSymmetric(result@aggregatedGuideGIs@data[,, replicate_name]))
  }
  expect_identical(dim(result@guideLFCs@data), dim(result@guideGIs@data))
  expect_identical(
    dimnames(result@guideLFCs@data),
    dimnames(result@guideGIs@data)
  )
  expected <- mean(c(
    input$GI[input$query_gene == "G1" & input$library_gene == "G2"],
    input$GI[input$query_gene == "G2" & input$library_gene == "G1"]
  ))
  expect_equal(
    result@guideGIs@data["G1", "G2", 1],
    input$GI[input$query_gene == "G1" & input$library_gene == "G2"]
  )
  expect_equal(result@aggregatedGuideGIs@data["G1", "G2", 1], expected)

  retained <- GIScores(
    input,
    pos_agnostic = TRUE,
    retain_directional = TRUE,
    block_layer = "guide_pair"
  )
  expect_true(retained@metadata$retain_directional)
})

test_that("GIScores creates a global position-agnostic pair matrix", {
  input <- make_multiplex_scores()
  input$GI <- as.numeric(factor(input$query_gene)) *
    100 +
    as.numeric(factor(input$library_gene))
  input$LFC <- input$GI / 10

  result <- GIScores(
    input,
    pos_agnostic = TRUE,
    symmetric_analysis_method = "global_preaverage",
    block_layer = "guide_pair"
  )

  expect_s4_class(result, "PosAgnMultiplexScreen")
  expect_identical(
    result@metadata$symmetric_analysis_method,
    "global_preaverage"
  )
  expect_identical(result@guideGIs@space, c("query_gene", "library_gene"))
  expect_identical(result@aggregatedGuideGIs@space, "gene_pair")
  expect_equal(
    rownames(result@aggregatedGuideGIs@data),
    result@screen_attr$unique_pairs
  )
  expect_equal(
    ncol(result@aggregatedGuideGIs@data),
    length(result@aggregatedGuideGIs@block_description)
  )
  expect_equal(
    length(result@aggregatedGuideGIs@blocks),
    ncol(result@aggregatedGuideGIs@data)
  )

  expected <- mean(c(
    input$GI[input$query_gene == "G1" & input$library_gene == "G2"],
    input$GI[input$query_gene == "G2" & input$library_gene == "G1"]
  ))

  expect_equal(result@aggregatedGuideGIs@data["G1;G2", 1L], expected)
  expect_identical(result@guideLFCs@space, result@guideGIs@space)
  expect_identical(dim(result@guideLFCs@data), dim(result@guideGIs@data))
  expect_identical(
    dimnames(result@guideLFCs@data),
    dimnames(result@guideGIs@data)
  )
})
