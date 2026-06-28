setMethod(
  "screenReport",
  signature = signature(GI_obj = "ScreenBase"),
  function(GI_obj, interactive = FALSE, print = TRUE, width = 80) {
    .line <- function(char = "-") paste(rep(char, width), collapse = "")

    .value <- function(x, default = "not available") {
      if (is.null(x) || length(x) == 0 || all(is.na(x))) {
        return(default)
      }
      paste(x, collapse = ", ")
    }

    .flag <- function(x) {
      if (isTRUE(x)) {
        "OK"
      } else if (identical(x, FALSE)) {
        "PROBLEM"
      } else {
        "not checked"
      }
    }

    .show_n <- function(x, n = 8) {
      if (is.null(x) || length(x) == 0) {
        return("none")
      }
      x <- as.character(x)
      suffix <- if (length(x) > n) {
        paste0(" ... +", length(x) - n, " more")
      } else {
        ""
      }
      paste0(paste(utils::head(x, n), collapse = ", "), suffix)
    }

    .class <- class(GI_obj)[1]
    .type <- c(
      "FixedPairScreen" = "fixed-pair",
      "MultiplexScreen" = "multiplex",
      "PosAgnMultiplexScreen" = "position-agnostic multiplex"
    )
    .type <- .value(.type[.class], .class)

    .a <- GI_obj@screen_attr
    .checks <- GI_obj@checks
    .md <- GI_obj@metadata
    .err <- GI_obj@errors
    .failed_queries <- .err$query_genes_not_usable
    if (is.null(.failed_queries)) {
      .failed_queries <- character()
    }

    .overview <- c(
      paste0("Screen class: ", .class),
      paste0("Interpreted design: ", .type),
      paste0("Query genes: ", .value(.a$n_query_genes)),
      paste0("Library genes: ", .value(.a$n_lib_genes)),
      paste0("All genes: ", .value(.a$n_all_genes)),
      paste0("Observed directional pairs: ", .value(length(.a$all_pairs))),
      paste0("Observed unordered pairs: ", .value(length(.a$unique_pairs))),
      paste0(
        "Guide-level replicate layers: ",
        .value(GI_obj@guideGIs@replicates)
      ),
      paste0(
        "Duplicate-correlation block layer: ",
        .value(GI_obj@guideGIs@block_layer)
      )
    )

    .decisions <- c(
      paste0("Selected model strategy: ", .type),
      paste0("Forced fixed-pair run: ", isTRUE(.md$force_fixed_pair)),
      paste0(
        "Position-agnostic output: ",
        methods::is(GI_obj, "PosAgnMultiplexScreen")
      ),
      paste0("Model objects available: ", length(GI_obj@limma_models)),
      paste0("Gene-level GI array available: ", length(GI_obj@geneGIs) > 0)
    )

    .check_lines <- if (length(.checks) == 0) {
      "No screen checks stored."
    } else {
      c(
        paste0(
          "Gene sets sufficiently overlapping: ",
          .flag(.checks$gene_sets_equal)
        ),
        paste0("Query size sufficient: ", .flag(.checks$query_sufficient)),
        paste0("Library size sufficient: ", .flag(.checks$library_sufficient)),
        paste0(
          "Stable library size across queries: ",
          .flag(.checks$stable_library_size)
        ),
        paste0(
          "Sufficient tests per query: ",
          .flag(.checks$sufficient_tests_per_query)
        ),
        paste0("Median tests per query: ", .value(.checks$avg_tests_per_query))
      )
    }

    .model_errors <- purrr::compact(.err$GI_computation_errors)
    .problem_lines <- c(
      paste0(
        "Query genes missing from library: ",
        length(.a$query_genes_not_in_lib),
        if (length(.a$query_genes_not_in_lib) > 0) {
          paste0(" (", .show_n(.a$query_genes_not_in_lib), ")")
        } else {
          ""
        }
      ),
      paste0(
        "Library genes missing from query set: ",
        length(.a$library_genes_not_in_query),
        if (length(.a$library_genes_not_in_query) > 0) {
          paste0(" (", .show_n(.a$library_genes_not_in_query), ")")
        } else {
          ""
        }
      ),
      paste0(
        "Query genes without usable model: ",
        length(.failed_queries)
      ),
      paste0("Stored model errors: ", length(.model_errors))
    )

    if (length(.failed_queries) > 0) {
      .problem_lines <- c(
        .problem_lines,
        paste0("Affected query genes: ", .show_n(.failed_queries))
      )
    }

    .report <- list(
      overview = .overview,
      decisions = .decisions,
      checks = .check_lines,
      problems = .problem_lines
    )

    if (isTRUE(print)) {
      cat("\n", .line("="), "\n", sep = "")
      cat("CeRberus screen report\n")
      cat(.line("="), "\n", sep = "")

      purrr::iwalk(.report, function(.x, .y) {
        cat("\n", toupper(.y), "\n", sep = "")
        cat(.line("-"), "\n", sep = "")
        cat(paste0("- ", .x, collapse = "\n"), "\n", sep = "")
      })
    }

    if (isTRUE(interactive)) {
      return(.report)
    }

    invisible(.report)
  }
)
