read_instructions <- function(yaml_fpath) {
  stopifnot(
    "yaml_fpath must be a single path string." = is.character(yaml_fpath) &&
      length(yaml_fpath) == 1 &&
      !is.na(yaml_fpath) &&
      nzchar(yaml_fpath),
    "Instruction file does not exist!" = file.exists(yaml_fpath)
  )

  instr <- read_yaml(file = yaml_fpath)
  if (is.null(instr) || !is.list(instr)) {
    stop("Instruction file must contain a YAML mapping/object.", call. = FALSE)
  }

  as_bool <- function(x, default = FALSE) {
    if (is.null(x)) {
      return(default)
    }
    if (is.logical(x)) {
      return(isTRUE(x))
    }
    vals <- c(
      "TRUE" = TRUE,
      "T" = TRUE,
      "True" = TRUE,
      "yes" = TRUE,
      "FALSE" = FALSE,
      "F" = FALSE,
      "False" = FALSE,
      "no" = FALSE
    )
    if (!as.character(x) %in% names(vals)) {
      return(default)
    }
    unname(vals[as.character(x)])
  }

  stopifnot(
    "Please provide a scores file using the 'scores_file' argument!" = "scores_file" %in%
      names(instr),
    "scores_file must be a single path string." = is.character(
      instr$scores_file
    ) &&
      length(instr$scores_file) == 1 &&
      !is.na(instr$scores_file) &&
      nzchar(instr$scores_file),
    "Given scores file does not exist!" = file.exists(instr$scores_file),
    "Output directory required!" = "output_directory" %in% names(instr),
    "output_directory must be a single path string." = is.character(
      instr$output_directory
    ),
    "Output directory must not be empty!" = length(instr$output_directory) ==
      1 &&
      !is.na(instr$output_directory) &&
      nzchar(instr$output_directory)
  )

  if (
    !"FDR" %in% names(instr) ||
      !instr$FDR %in% c("BH", "bonferroni")
  ) {
    instr$FDR <- "BH"
  }

  #  if (!"output_prefix" %in% names(instr) ||
  #      instr$output_prefix == "") {
  #    instr$output_prefix <- gsub(".csv$|.rds$", "",
  #                                basename(instr$scores_file))
  #  }

  instr$overwrite_output <- as_bool(instr$overwrite_output, default = TRUE)
  instr$verbose <- as_bool(instr$verbose, default = FALSE)

  instr$output_directory <- normalizePath(
    instr$output_directory,
    winslash = "/",
    mustWork = FALSE
  )

  if (grepl("[/\\\\]$", instr$output_directory)) {
    warning("Removing trailing path separator from output_directory.")
    instr$output_directory <- sub("[/\\\\]+$", "", instr$output_directory)
  }

  return(instr)
}


collect_all_layer_configurations <- function(
  GI_data,
  .to_use = c("tech_rep", "bio_rep", "guide_pair"),
  make_pos_agnostic,
  verbose = FALSE
) {
  .collapsable_layers <- intersect(
    c("tech_rep", "bio_rep", "guide_pair"),
    colnames(GI_data)
  )
  .to_use <- intersect(.to_use, colnames(GI_data))

  output <- list()

  for (.use in .to_use) {
    output[[str_c("default_", .use, "_used")]] <- GIScores(
      GI_data,
      block_layer = .use,
      pos_agnostic = make_pos_agnostic,
      verbose = verbose
    )
  }

  for (.i in seq_len(length(.collapsable_layers) - 1)) {
    for (.to_collapse in combn(.collapsable_layers, .i, simplify = FALSE)) {
      for (.use in setdiff(.collapsable_layers, .to_collapse)) {
        .n <- str_c(
          str_c(.to_collapse, collapse = "_"),
          "_collapsed_",
          .use,
          "_used"
        )
        if (isTRUE(verbose)) {
          print(.n)
        }
        output[[.n]] <- GIScores(
          GI_data,
          collapse_layers = .to_collapse,
          block_layer = .use,
          pos_agnostic = make_pos_agnostic,
          verbose = verbose
        )
      }
    }
  }
  return(output)
}


find_optimal_configuration <- function(GI_list, keep_all = FALSE) {
  stopifnot(
    "GI_list must contain at least one screen object." = length(GI_list) > 0,
    "keep_all must be TRUE or FALSE." = is.logical(keep_all) &&
      length(keep_all) == 1 &&
      !is.na(keep_all)
  )

  .x <- GI_list |> map(~ .x@dupCorrelation) |> map_dbl(mean, na.rm = TRUE)

  .d <- data.table(
    config = names(GI_list),
    dcor = GI_list |>
      map(~ .x@dupCorrelation) |>
      map_dbl(mean, na.rm = TRUE) |>
      round(3)
  )

  if (any(.x >= 0)) {
    .x <- keep(.x, .x >= 0)
  } # if greater or = 0 found, remove <= 0
  if (any(.x > 0)) {
    .x <- .x |> keep(.x > 0)
  } # if any greater 0 found, remove = 0
  if (length(.x) >= 1) {
    .x <- .x[which.min(.x)]
  } # if more than one option remains, choose weakest correlation

  .d[,
    kept := fcase(
      config %in% names(.x) , "selected" ,
      default = ""
    )
  ]

  for (.n in names(GI_list)) {
    GI_list[[.n]]@metadata$dupcor_data <- .d
  }

  if (keep_all) {
    return(GI_list)
  }
  if (!keep_all) {
    return(GI_list[names(.x)])
  }
}
