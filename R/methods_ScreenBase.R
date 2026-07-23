# gGRNA_GIs methods

#purrr::walk(c("guideGIs"), ~ {
#  .x2 <- paste0(.x, "<-")
#  setMethod(.x, "guideGIs", .f1))
#  setMethod(.x2, "guideGIs", .f2)
#}
#)

# ScreenBase methods

setMethod("checks", "ScreenBase", function(x) {
  return(slot(x, "checks"))
})

#####

setMethod("checks<-", "ScreenBase", function(x, value) {
  slot(x, "checks") <- value
  return(x)
})

#####

setMethod("dupCorrelation", "ScreenBase", function(x) {
  return(slot(x, "dupCorrelation"))
})

#####

setMethod("dupCorrelation<-", "ScreenBase", function(x, value) {
  slot(x, "dupCorrelation") <- value
  return(x)
})

#####

setMethod("errors", "ScreenBase", function(x) {
  return(slot(x, "errors"))
})

#####

setMethod("errors<-", "ScreenBase", function(x, value) {
  slot(x, "errors") <- value
  return(x)
})

#####

setMethod("geneGIs", "ScreenBase", function(x) {
  return(slot(x, "geneGIs"))
})

#####

setMethod("geneGIs<-", "ScreenBase", function(x, value) {
  slot(x, "geneGIs") <- value
  return(x)
})

#####

setMethod("guideGIs", "ScreenBase", function(x) {
  return(slot(x, "guideGIs"))
})

#####

setMethod("guideGIs<-", "ScreenBase", function(x, value) {
  slot(x, "guideGIs") <- value
  return(x)
})

#####

setMethod("limma_models", "ScreenBase", function(x) {
  return(slot(x, "limma_models"))
})

#####

setMethod("limma_models<-", "ScreenBase", function(x, value) {
  slot(x, "limma_models") <- value
  return(x)
})

#####

setMethod("screen_attr", "ScreenBase", function(x) {
  return(slot(x, "screen_attr"))
})

#####

setMethod("screen_attr<-", "ScreenBase", function(x, value) {
  slot(x, "screen_attr") <- value
  return(x)
})

#####

setMethod("symmGeneGIs", "ScreenBase", function(x) {
  stop(
    "symmGeneGIs is only available for PosAgnMultiplexScreen objects. ",
    "Run with pos_agnostic = TRUE first.",
    call. = FALSE
  )
})

#####

setMethod("symmGeneGIs<-", "ScreenBase", function(x, value) {
  stop(
    "symmGeneGIs can only be assigned for PosAgnMultiplexScreen objects.",
    call. = FALSE
  )
})

#####

setMethod("symmGeneGIs", "PosAgnMultiplexScreen", function(x) {
  slot(x, "symmGeneGIs")
})

#####

setMethod("symmGeneGIs<-", "PosAgnMultiplexScreen", function(x, value) {
  slot(x, "symmGeneGIs") <- value
  x
})

#####

setMethod(
  "import_scores",
  signature = signature(GI_obj = "ScreenBase"),
  function(GI_obj) {
    .md <- GI_obj@metadata

    .md$input <- copy(.md$input)

    stopifnot(
      "The input object needs to be a data frame or data table." = is.data.frame(
        .md$input
      ),
      "The query gene column is not in the provided dataset." = .md$query_col %in%
        colnames(.md$input),
      "The library gene column is not in the provided dataset." = .md$lib_col %in%
        colnames(.md$input)
    )

    .md$input <- as.data.table(.md$input)

    setnames(
      .md$input,
      old = c(
        .md$bio_rep_col,
        .md$tech_rep_col,
        .md$guide_col,
        .md$query_col,
        .md$lib_col,
        .md$gi_col
      ),
      new = c(
        "bio_rep",
        "tech_rep",
        "guide_pair",
        "query_gene",
        "library_gene",
        "GI"
      ),
      skip_absent = TRUE
    )

    .md$input[,
      gene_pair := paste0(get("query_gene"), ";", get("library_gene"))
    ]

    #            .md$input[, replicate := do.call(paste, c(.SD, sep = "_")),
    #                      .SDcols = intersect(c("bio_rep", "tech_rep", "guide_pair"), colnames(.md$input))]

    GI_obj@metadata <- .md

    return(GI_obj)
  }
)

#####

setMethod(
  "get_screen_attributes",
  signature = signature(GI_obj = "ScreenBase"),
  function(GI_obj) {
    .i <- GI_obj@metadata$input
    .a <- list(contrasts = NULL)
    .a$query_genes <- .i[, unique(get("query_gene"))]
    .a$library_genes <- .i[, unique(get("library_gene"))]
    .a$all_genes <- union(.a$query_genes, .a$library_genes)
    .a$query_genes_not_in_lib <- setdiff(.a$query_genes, .a$library_genes)
    .a$library_genes_not_in_query <- setdiff(.a$library_genes, .a$query_genes)

    .a$n_query_genes <- length(.a$query_genes)
    .a$n_lib_genes <- length(.a$library_genes)
    .a$n_all_genes <- length(.a$all_genes)

    .a$observations_per_query <- purrr::map_int(
      purrr::set_names(.a$query_genes),
      \(.g) {
        .i[query_gene == .g, .N]
      }
    )

    .a$all_pairs <- .i[, unique(get("gene_pair"))]
    .a$unique_pairs <- .i[, unique(sort_gene_pairs(
      get("query_gene"),
      get("library_gene")
    ))]

    GI_obj@screen_attr <- .a
    return(GI_obj)
  }
)

#####

setMethod(
  "run_checks",
  signature = signature(GI_obj = "ScreenBase"),
  function(GI_obj, min_query_size = 20, min_library_size = 20) {
    .a <- GI_obj@screen_attr

    GI_obj@checks <- list(
      gene_sets_equal = (length(.a$query_genes_not_in_lib) <=
        0.02 * .a$n_query_genes) &
        (length(.a$library_genes_not_in_query) <= 0.02 * .a$n_lib_genes),
      query_sufficient = .a$n_query_genes >= min_query_size,
      library_sufficient = .a$n_lib_genes >= min_library_size,
      stable_library_size = sum(
        .a$observations_per_query !=
          median(.a$observations_per_query, na.rm = TRUE)
      ) <=
        10,
      sufficient_tests_per_query = sum(
        .a$observations_per_query >= min_library_size
      ) >=
        0.95 * length(.a$observations_per_query),
      avg_tests_per_query = median(.a$observations_per_query, na.rm = TRUE)
    )
    return(GI_obj)
  }
)

#####

setMethod(
  "set_screenType",
  signature = signature(GI_obj = "ScreenBase"),
  function(GI_obj) {
    .md <- GI_obj@metadata
    .c <- GI_obj@checks

    .is_multiplex_candidate <- isTRUE(.c$query_sufficient) &&
      isTRUE(.c$library_sufficient) &&
      #isTRUE(.c$gene_sets_equal) &&
      isTRUE(.c$sufficient_tests_per_query)

    .inferred_screen_type <- if (.is_multiplex_candidate) {
      "multiplex"
    } else {
      "fixed_pair"
    }

    .requested_screen_type <- .md$requested_screen_type
    if (is.null(.requested_screen_type)) {
      .requested_screen_type <- "auto"
    }

    .selected_screen_type <- if (identical(.requested_screen_type, "auto")) {
      .inferred_screen_type
    } else {
      .requested_screen_type
    }

    if (.is_multiplex_candidate) {
      if (!isTRUE(.c$stable_library_size)) {
        warning(
          "Multiplex screen has variable observations per query; ",
          "keeping multiplex classification.",
          call. = FALSE
        )
      }
    }

    if (!identical(.selected_screen_type, .inferred_screen_type)) {
      warning(
        "Requested screen_type '",
        .selected_screen_type,
        "' overrides inferred screen type '",
        .inferred_screen_type,
        "'.",
        call. = FALSE
      )
    }

    .type <- switch(
      .selected_screen_type,
      fixed_pair = "FixedPairScreen",
      multiplex = "MultiplexScreen"
    )

    .md$requested_screen_type <- .requested_screen_type
    .md$inferred_screen_type <- .inferred_screen_type
    .md$selected_screen_type <- .selected_screen_type

    GI_obj <- as(object = GI_obj, Class = .type)

    if (is(GI_obj, "MultiplexScreen")) {
      GI_obj@guideGIs@space <- c("query_gene", "library_gene")
    } else if (is(GI_obj, "FixedPairScreen")) {
      GI_obj@guideGIs@space <- "gene_pair"
    }

    GI_obj@guideGIs@replicates <- intersect(
      c("guide_pair", "tech_rep", "bio_rep"),
      colnames(.md$input)
    )

    GI_obj@metadata <- .md
    return(GI_obj)
  }
)

#####

setMethod(
  "compute_dupCorrelation",
  signature = signature(.x = "ScreenBase"),
  function(.x) {
    .x@dupCorrelation <- compute_dupCorrelation(.x@guideGIs)
    return(.x)
  }
)

#####

setMethod(
  "compute_models",
  signature = signature(GI_obj = "FixedPairScreen"),
  function(GI_obj) {
    GI_obj@limma_models <- limma::lmFit(
      object = GI_obj@guideGIs@data,
      block = GI_obj@guideGIs@blocks,
      correlation = GI_obj@dupCorrelation
    ) |>
      limma::eBayes()
    return(GI_obj)
  }
)

#####

setMethod(
  "compute_models",
  signature = signature(GI_obj = "MultiplexScreen"),
  function(GI_obj) {
    output <- GI_obj@screen_attr$query_genes |>
      set_names()
    output <- output |>
      map(safely(\(.g) {
        .fit <- limma::lmFit(
          object = GI_obj@guideGIs@data[.g, GI_obj@screen_attr$library_genes, ],
          block = GI_obj@guideGIs@blocks,
          correlation = GI_obj@dupCorrelation[[.g]]
        ) |>
          limma::eBayes()

        return(.fit)
      }))

    GI_obj@limma_models <- map(output, "result")
    GI_obj@errors$query_genes_not_usable <- map(output, "result") |>
      keep(is.null) |>
      names()
    GI_obj@errors$GI_computation_errors <- map(output, "error")

    #  map("error") |>
    #  compact() |>
    #  map_chr(~ .x$message)

    if (length(GI_obj@errors$query_genes_not_usable) > 0) {
      warning(str_c(
        "Failed computing GIs for ",
        length(GI_obj@errors$query_genes_not_usable),
        " genes."
      ))
    }
    return(GI_obj)
  }
)

#####

setMethod(
  "compute_models",
  signature = signature(GI_obj = "PosAgnMultiplexScreen"),
  function(GI_obj) {
    symmetric_analysis_method <- get_symmetric_analysis_method(GI_obj)

    if (identical(symmetric_analysis_method, "preaverage")) {
      return(methods::callNextMethod(GI_obj))
    }

    if (!identical(symmetric_analysis_method, "global_preaverage")) {
      stop(
        "Unsupported symmetric analysis method: ",
        symmetric_analysis_method,
        call. = FALSE
      )
    }

    if (length(dim(GI_obj@guideGIs@data)) != 2L) {
      stop(
        "global_preaverage requires a pair-by-observation guide-level matrix.",
        call. = FALSE
      )
    }

    if (
      !identical(
        rownames(GI_obj@guideGIs@data),
        GI_obj@screen_attr$unique_pairs
      )
    ) {
      stop(
        "The global model matrix rows do not match the stored unordered pairs.",
        call. = FALSE
      )
    }

    if (length(GI_obj@dupCorrelation) != 1L) {
      stop(
        "global_preaverage requires one global duplicate-correlation estimate.",
        call. = FALSE
      )
    }

    if (
      isTRUE(GI_obj@guideGIs@use_blocks) &&
        length(GI_obj@guideGIs@blocks) != ncol(GI_obj@guideGIs@data)
    ) {
      stop(
        "The number of block assignments does not match the global model columns.",
        call. = FALSE
      )
    }

    GI_obj@limma_models <- limma::lmFit(
      object = GI_obj@guideGIs@data,
      block = GI_obj@guideGIs@blocks,
      correlation = GI_obj@dupCorrelation
    ) |>
      limma::eBayes()

    return(GI_obj)
  }
)

#####

setMethod(
  "collect_GIs",
  signature = signature(GI_obj = "FixedPairScreen"),
  function(GI_obj, FDR_method = "BH") {
    stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)

    output <- list(rownames(GI_obj@guideGIs@data), c("GI", "pval", "FDR"))

    output <- array(
      data = NA,
      dim = purrr::map_int(output, length),
      dimnames = output
    )

    output[, "GI"] <- GI_obj@limma_models$coefficients[, 1]
    output[, "pval"] <- GI_obj@limma_models$p.value[, 1]
    output[, "FDR"] <- stats::p.adjust(
      GI_obj@limma_models$p.value[, 1],
      method = FDR_method
    )

    GI_obj@geneGIs <- output

    return(GI_obj)
  }
)

#####

setMethod(
  "collect_GIs",
  signature = signature(GI_obj = "MultiplexScreen"),
  function(GI_obj, FDR_method = "BH") {
    stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)

    output <- GI_obj@limma_models |>
      imap(\(.m, .y) {
        .x <- data.frame(
          library_gene = GI_obj@screen_attr$library_genes,
          query_gene = .y
        )
        if (is.null(.m)) {
          .x$GI <- NA
          .x$pval <- NA
          .x$FDR <- NA
        } else {
          .x$GI <- .m$coefficients[, 1]
          .x$pval <- .m$p.value[, 1]
          .x$FDR <- stats::p.adjust(.m$p.value[, 1], method = FDR_method)
        }
        return(.x)
      })

    output <- output |>
      data.table::rbindlist(fill = TRUE) |>
      data.table::melt.data.table(measure.vars = c("GI", "pval", "FDR")) |>
      reshape2::acast(
        formula = as.formula("query_gene ~ library_gene ~ variable"),
        value.var = "value",
        drop = FALSE
      )

    GI_obj@geneGIs <- output

    if (
      length(setdiff(
        rownames(GI_obj@guideGIs@data),
        rownames(GI_obj@geneGIs)
      )) !=
        0 |
        length(setdiff(
          colnames(GI_obj@guideGIs@data),
          colnames(GI_obj@geneGIs)
        )) !=
          0
    ) {
      warning("Some genes were lost!")
    }
    return(GI_obj)
  }
)

#####

setMethod(
  "collect_GIs",
  signature = signature(GI_obj = "PosAgnMultiplexScreen"),
  function(GI_obj, FDR_method = "BH") {
    stopifnot(
      "Unknown FDR method provided." = FDR_method %in% p.adjust.methods
    )

    symmetric_analysis_method <- get_symmetric_analysis_method(GI_obj)

    if (identical(symmetric_analysis_method, "preaverage")) {
      GI_obj <- methods::callNextMethod(GI_obj)

      .x <- data.table(gene_pair = GI_obj@screen_attr$unique_pairs)

      .x[, query_gene := str_split_i(gene_pair, ";", 1)]
      .x[, library_gene := str_split_i(gene_pair, ";", 2)]
      .x[,
        GI := gather_symmetric_scores(
          pairs = gene_pair,
          .arr = GI_obj@geneGIs[,, "GI"]
        )
      ]
      .x[, GI_z := z_transform(GI)]
      .x[,
        pval := gather_symmetric_scores(
          pairs = gene_pair,
          .arr = GI_obj@geneGIs[,, "pval"]
        )
      ]
      .x[,
        FDR := balanced_FDR(
          pairs = gene_pair,
          pval_array = GI_obj@geneGIs[,, "pval"],
          fdr_method = FDR_method
        )
      ]

      GI_obj@symmGeneGIs <- .x

      return(GI_obj)
    }

    if (!identical(symmetric_analysis_method, "global_preaverage")) {
      stop(
        "Unsupported symmetric analysis method: ",
        symmetric_analysis_method,
        call. = FALSE
      )
    }

    model_pairs <- rownames(GI_obj@limma_models$coefficients)
    expected_pairs <- GI_obj@screen_attr$unique_pairs

    if (!identical(model_pairs, expected_pairs)) {
      stop(
        "The global limma output rows do not match the stored unordered pairs.",
        call. = FALSE
      )
    }

    GI <- GI_obj@limma_models$coefficients[, 1L]
    pval <- GI_obj@limma_models$p.value[, 1L]
    FDR <- stats::p.adjust(pval, method = FDR_method)

    GI_obj@geneGIs <- cbind(
      GI = as.numeric(GI),
      pval = as.numeric(pval),
      FDR = as.numeric(FDR)
    )

    rownames(GI_obj@geneGIs) <- model_pairs

    GI_obj@symmGeneGIs <- data.table::data.table(
      gene_pair = model_pairs,
      query_gene = stringr::str_split_i(model_pairs, ";", 1),
      library_gene = stringr::str_split_i(model_pairs, ";", 2),
      GI = as.numeric(GI),
      GI_z = z_transform(as.numeric(GI)),
      pval = as.numeric(pval),
      FDR = as.numeric(FDR)
    )

    return(GI_obj)
  }
)

#####

setMethod(
  "GI_df",
  signature = signature(GI_obj = "FixedPairScreen"),
  function(GI_obj) {
    output <- GI_obj@geneGIs

    output <- data.table::data.table(
      gene_pair = rownames(output),
      query_gene = str_split_i(rownames(output), ";", 1),
      library_gene = str_split_i(rownames(output), ";", 2),
      output
    )
    return(output)
  }
)

#####

setMethod(
  "GI_df",
  signature = signature(GI_obj = "MultiplexScreen"),
  function(GI_obj) {
    output <- GI_obj@geneGIs

    output <- output |>
      flatten_array(
        dnames = c("query_gene", "library_gene", "variable"),
        value.name = "value"
      ) |>
      data.table::dcast(
        formula = query_gene + library_gene ~ variable,
        value.var = "value"
      )
    output[, gene_pair := str_c(query_gene, ";", library_gene)]
    output <- output[,
      .SD,
      .SDcols = c(
        "gene_pair",
        "query_gene",
        "library_gene",
        "GI",
        "pval",
        "FDR"
      )
    ]

    return(output)
  }
)

#####

setMethod(
  "GI_df",
  signature = signature(GI_obj = "PosAgnMultiplexScreen"),
  function(GI_obj) {
    GI_obj@symmGeneGIs
  }
)

#####

setMethod(
  "dupCorrelation_df",
  signature = signature(GI_obj = "ScreenBase"),
  function(GI_obj) {
    output <- data.table(data.frame(dupcor = GI_obj@dupCorrelation))
    return(output)
  }
)

#####

setMethod(
  "create_log",
  signature = signature(GI_obj = "ScreenBase"),
  function(
    GI_obj,
    status = "available",
    stage = "not specified",
    condition = NULL,
    max_items = 10L
  ) {
    stopifnot(
      "status must be a single non-missing character value." = is.character(
        status
      ) &&
        length(status) == 1L &&
        !is.na(status),
      "stage must be a single non-missing character value." = is.character(
        stage
      ) &&
        length(stage) == 1L &&
        !is.na(stage),
      "max_items must be a single positive whole number." = is.numeric(
        max_items
      ) &&
        length(max_items) == 1L &&
        !is.na(max_items) &&
        is.finite(max_items) &&
        max_items >= 1L &&
        max_items == as.integer(max_items)
    )

    max_items <- as.integer(max_items)

    value_or <- function(x, default = "not available") {
      if (is.null(x) || length(x) == 0L || all(is.na(x))) {
        return(default)
      }
      paste(as.character(x), collapse = ", ")
    }

    show_items <- function(x) {
      if (is.null(x) || length(x) == 0L) {
        return("none")
      }
      x <- as.character(x)
      shown <- utils::head(x, max_items)
      suffix <- if (length(x) > max_items) {
        paste0(" ... +", length(x) - max_items, " more")
      } else {
        ""
      }
      paste0(paste(shown, collapse = " | "), suffix)
    }

    array_summary <- function(x) {
      if (length(x) == 0L) {
        return("empty")
      }
      dims <- dim(x)
      dim_text <- if (is.null(dims)) {
        length(x)
      } else {
        paste(dims, collapse = " x ")
      }
      paste0(
        "dimensions=",
        dim_text,
        "; values=",
        length(x),
        "; missing=",
        sum(is.na(x))
      )
    }

    collect_error_messages <- function(x) {
      if (is.null(x)) {
        return(character())
      }
      if (inherits(x, "condition")) {
        return(conditionMessage(x))
      }
      if (is.list(x)) {
        return(unlist(lapply(x, collect_error_messages), use.names = FALSE))
      }
      if (is.character(x)) {
        return(x[!is.na(x) & nzchar(x)])
      }
      character()
    }

    condition_message <- if (is.null(condition)) {
      "none"
    } else if (inherits(condition, "condition")) {
      conditionMessage(condition)
    } else {
      paste(as.character(condition), collapse = " | ")
    }

    package_version <- tryCatch(
      as.character(utils::packageVersion("CeRberus")),
      error = function(e) "development version"
    )

    .a <- GI_obj@screen_attr
    .checks <- GI_obj@checks
    .md <- GI_obj@metadata
    .errors <- GI_obj@errors
    .guide <- GI_obj@guideGIs

    dupcor <- GI_obj@dupCorrelation
    finite_dupcor <- dupcor[is.finite(dupcor)]
    dupcor_summary <- if (length(dupcor) == 0L) {
      "not available"
    } else {
      paste0(
        "values=",
        length(dupcor),
        "; finite=",
        length(finite_dupcor),
        "; missing=",
        sum(is.na(dupcor)),
        if (length(finite_dupcor) > 0L) {
          paste0(
            "; min=",
            signif(min(finite_dupcor), 4L),
            "; median=",
            signif(stats::median(finite_dupcor), 4L),
            "; max=",
            signif(max(finite_dupcor), 4L)
          )
        } else {
          ""
        }
      )
    }

    model_count <- if (inherits(GI_obj@limma_models, "MArrayLM")) {
      1L
    } else {
      sum(!purrr::map_lgl(GI_obj@limma_models, is.null))
    }

    check_lines <- if (length(.checks) == 0L) {
      "- none recorded"
    } else {
      paste0(
        "- ",
        names(.checks),
        "=",
        purrr::map_chr(.checks, value_or),
        collapse = "\n"
      )
    }

    failed_queries <- .errors$query_genes_not_usable
    stored_errors <- collect_error_messages(.errors$GI_computation_errors)

    lines <- c(
      "CeRberus diagnostic log",
      "========================",
      paste0("Timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      paste0("Status: ", status),
      paste0("Pipeline stage: ", stage),
      paste0("Condition: ", condition_message),
      "",
      "Runtime",
      "-------",
      paste0("CeRberus version: ", package_version),
      paste0("R version: ", R.version.string),
      paste0("Platform: ", R.version$platform),
      paste0("Operating system: ", value_or(Sys.info()["sysname"])),
      "",
      "Screen",
      "------",
      paste0("Configuration: ", value_or(.md$config)),
      paste0("Class: ", class(GI_obj)[1L]),
      paste0("Query genes: ", value_or(.a$n_query_genes)),
      paste0("Library genes: ", value_or(.a$n_lib_genes)),
      paste0("All genes: ", value_or(.a$n_all_genes)),
      paste0("Directional pairs: ", value_or(length(.a$all_pairs))),
      paste0("Unordered pairs: ", value_or(length(.a$unique_pairs))),
      paste0("Guide GI data: ", array_summary(.guide@data)),
      paste0("Gene GI data: ", array_summary(GI_obj@geneGIs)),
      "",
      "Replicates and modelling",
      "------------------------",
      paste0("Replicate layers: ", value_or(.guide@replicates)),
      paste0("Collapsed layers: ", value_or(.guide@collapse, "none")),
      paste0("Block layer: ", value_or(.guide@block_layer, "none")),
      paste0("Uses blocks: ", value_or(.guide@use_blocks)),
      paste0("Block assignments: ", length(.guide@blocks)),
      paste0("Duplicate correlation: ", dupcor_summary),
      paste0("Fitted models: ", model_count),
      "",
      "Checks",
      "------",
      check_lines,
      "",
      "Errors and incomplete results",
      "-----------------------------",
      paste0(
        "Failed queries (",
        length(failed_queries),
        "): ",
        show_items(failed_queries)
      ),
      paste0(
        "Stored model errors (",
        length(stored_errors),
        "): ",
        show_items(stored_errors)
      ),
      paste0("Gene-level results available: ", length(GI_obj@geneGIs) > 0L)
    )

    paste(lines, collapse = "\n")
  }
)

#####

setMethod(
  "symmetry_test",
  signature = signature(GI_obj = "MultiplexScreen"),
  function(GI_obj, cutoff = 0.99) {
    .test <- map_lgl(set_names(GI_obj@guideGIs@replicates), \(.r) {
      .x <- GI_obj@guideGIs@data[,, .r]
      if (all(.x == t(.x)) || all(dplyr::near(.x, t(.x)))) {
        return(TRUE)
      } else {
        return(
          cor(as.vector(.x), as.vector(t(.x)), use = "pairwise.complete.obs") >=
            cutoff
        )
      }
    })

    return(all(.test))
  }
)

#####
