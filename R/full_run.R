#' Run the Full GI Analysis Pipeline.
#'
#' @param yaml_fpath Path to a YAML instruction file.
#' @param return_output Binary that defines if the function should save
#' @returns Saves the GI scores and plots to the output directory defined in
#' the YAML file.
#' @export

full_run <- function(yaml_fpath, return_output = TRUE) {
  instr <- read_instructions(yaml_fpath)

  stopifnot("Output directory required!" = !is.null(instr$output_directory))

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

  .make_symmetric <- "make_symmetric" %in%
    names(instr) &&
    (isTRUE(instr$make_symmetric))

  .data <- collect_all_layer_configurations(
    .data,
    make_pos_agnostic = .make_symmetric,
    verbose = instr$verbose
  )

  .data <- map(.data, compute_dupCorrelation)

  # store all GI objects
  if (
    "output_directory" %in% names(instr) && (isTRUE(instr$overwrite_output))
  ) {
    saveRDS(.data, file.path(instr$output_directory, "all_GI_objects.rds"))
  }

  .keep_all <- "keep_all_configurations" %in%
    names(instr) &&
    (isTRUE(instr$keep_all_configurations))

  .data <- find_optimal_configuration(.data, keep_all = .keep_all)

  .data <- .data |>
    compute_dupcor_plot(
      .fpath = file.path(
        instr$output_directory,
        "duplicateCorrelationPlot.png"
      ),
      verbose = instr$verbose
    )

  .data <- map(.data, compute_models)

  .data <- map(.data, collect_GIs, FDR_method = instr$FDR)

  if (isTRUE(instr$verbose)) {
    purrr::iwalk(.data, function(.x, .y) {
      cat("\nConfiguration: ", .y, "\n", sep = "")
      screenReport(.x, interactive = FALSE, print = TRUE)
    })
  }

  if ("output_directory" %in% names(instr) & instr$overwrite_output) {
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

  #}####

  #if (instr$verbose) {print("(5/5)")}
  if (!return_output) {
    .data <- NULL
  }
  return(.data)
}
