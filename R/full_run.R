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

  .stage <- "read_scores"
  .data <- NULL

  write_pipeline_log <- function(
    path,
    status,
    stage,
    condition = NULL,
    screen_objects = NULL
  ) {
    tryCatch(
      {
        if (methods::is(screen_objects, "ScreenBase")) {
          screen_objects <- list(screen_objects)
        }

        if (!is.list(screen_objects)) {
          screen_objects <- list()
        }

        screen_objects <- screen_objects[
          purrr::map_lgl(
            screen_objects,
            ~ methods::is(.x, "ScreenBase")
          )
        ]

        if (length(screen_objects) == 0L) {
          screen_objects <- list(methods::new("ScreenBase"))
          names(screen_objects) <- "no screen object available"
        } else if (is.null(names(screen_objects))) {
          names(screen_objects) <- paste0(
            "configuration_",
            seq_along(screen_objects)
          )
        }

        log_blocks <- purrr::imap(
          screen_objects,
          function(.x, .name) {
            c(
              paste0("Configuration entry: ", .name),
              create_log(
                .x,
                status = status,
                stage = stage,
                condition = condition
              )
            )
          }
        )

        writeLines(
          unlist(
            Map(
              c,
              log_blocks,
              MoreArgs = list(c("", ""))
            ),
            use.names = FALSE
          ),
          con = path,
          useBytes = TRUE
        )
      },
      error = function(log_error) {
        warning(
          "Failed to write CeRberus diagnostic log: ",
          conditionMessage(log_error),
          call. = FALSE
        )
      }
    )
  }

  tryCatch(
    {
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

      .stage <- "construct_configurations"
      .data <- collect_all_layer_configurations(
        .data,
        screen_type = instr$screen_type,
        pos_agnostic = instr$pos_agnostic,
        symmetric_analysis_method = instr$symmetric_analysis_method,
        verbose = instr$verbose
      )

      .stage <- "duplicate_correlation"
      .data <- map(.data, compute_dupCorrelation)

      if (isTRUE(instr$overwrite_output)) {
        saveRDS(.data, file.path(instr$output_directory, "all_GI_objects.rds"))
      }

      .stage <- "select_configuration"
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

      .stage <- "fit_models"
      .data <- map(.data, compute_models)

      .stage <- "collect_GIs"
      .data <- map(.data, collect_GIs, FDR_method = instr$FDR)

      if (isTRUE(instr$verbose)) {
        purrr::iwalk(.data, function(.x, .y) {
          cat("\nConfiguration: ", .y, "\n", sep = "")
          screenReport(.x, interactive = FALSE, print = TRUE)
        })
      }

      .stage <- "write_outputs"
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

        write_pipeline_log(
          path = file.path(instr$output_directory, "CeRberus.log"),
          status = "completed",
          stage = "complete",
          screen_objects = .data
        )
      }

      if (!return_output) {
        return(NULL)
      }

      return(.data)
    },
    error = function(pipeline_error) {
      write_pipeline_log(
        path = file.path(instr$output_directory, "CeRberus_error.log"),
        status = "failed",
        stage = .stage,
        condition = pipeline_error,
        screen_objects = .data
      )
      stop(pipeline_error)
    }
  )
}

#####
