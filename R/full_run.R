#' Run the full GI analysis pipeline
#'
#' @param yaml_fpath Path to a YAML instruction file.
#' @param return_output Logical scalar. If `TRUE`, return the retained screen
#'   configurations. If `FALSE`, return `NULL` after completing the pipeline.
#'
#' @return A named list of retained CeRberus screen objects when
#'   `return_output = TRUE`; otherwise `NULL`.
#'
#' @export

#####

full_run <- function(yaml_fpath, return_output = TRUE) {
  stopifnot(
    "return_output must be TRUE or FALSE." = is.logical(return_output) &&
      length(return_output) == 1L &&
      !is.na(return_output)
  )

  instr <- read_instructions(yaml_fpath)

  dir.create(instr$output_directory, showWarnings = FALSE, recursive = TRUE)

  #for (.pkg in c("BiocManager")) {if (!require(.pkg, quietly = TRUE)) {utils::install.packages(.pkg)}}
  #if (!require("limma", quietly = TRUE)) {BiocManager::install("limma")}

  #if (instr$verbose) {print(str(instr, give.attr = FALSE))}

  scores_ext <- tolower(tools::file_ext(instr$scores_file))

  if (!nzchar(scores_ext)) {
    stop(
      "scores_file must have a file extension. Supported formats are: .csv and .rds.",
      call. = FALSE
    )
  }

  .data <- switch(
    scores_ext,
    "csv" = data.table::fread(instr$scores_file),
    "rds" = readRDS(instr$scores_file),
    stop(
      "Unsupported scores_file extension: .",
      scores_ext,
      ". Supported formats are: .csv and .rds.",
      call. = FALSE
    )
  )

  .data <- collect_all_layer_configurations(
    .data,
    pos_agnostic = instr$pos_agnostic,
    symmetric_analysis_method = instr$symmetric_analysis_method,
    verbose = instr$verbose
  )

  .data <- map(.data, compute_dupCorrelation)

  # store all GI objects
  if (
    "output_directory" %in% names(instr) && (isTRUE(instr$overwrite_output))
  ) {
    saveRDS(.data, file.path(instr$output_directory, "all_GI_objects.rds"))
  }

  .data <- find_optimal_configuration(
    .data,
    keep_all = instr$keep_all_configurations
  )

  if (instr$overwrite_output) {
    .data <- .data |>
      compute_dupcor_plot(
        .fpath = file.path(
          instr$output_directory,
          "duplicateCorrelationPlot.png"
        ),
        verbose = instr$verbose
      )
  }

  .data <- map(.data, compute_models)

  .data <- map(.data, collect_GIs, FDR_method = instr$FDR)

  if (isTRUE(instr$verbose)) {
    purrr::iwalk(.data, function(.x, .y) {
      cat("\nConfiguration: ", .y, "\n", sep = "")
      screenReport(.x, interactive = FALSE, print = TRUE)
    })
  }

  if (instr$overwrite_output) {
    .output <- list()

    .output$duplicate_correlation <- .data[[1]]@metadata$dupcor_data

    for (.n in names(.data)) {
      .output[[paste0("GI_scores_", .n)]] <- GI_df(.data[[.n]])
    }

    .output |>
      iwalk(
        ~ {
          data.table::fwrite(
            x = .x,
            file = file.path(instr$output_directory, paste0(.y, ".csv"))
          )
        }
      )

    #.log <- create_log(.data)

    #writeLines(.log, file.path(instr$output_directory, "limma.log"))
  }

  if (!return_output) {
    return(NULL)
  }

  return(.data)
}

#####
