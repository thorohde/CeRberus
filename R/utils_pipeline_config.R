#####

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

  if (!"symmetric_analysis_method" %in% names(instr)) {
    instr$symmetric_analysis_method <- "preaverage"
  }

  instr$symmetric_analysis_method <- validate_symmetric_analysis_method(
    instr$symmetric_analysis_method
  )

  instr$overwrite_output <- as_bool(instr$overwrite_output, default = TRUE)
  instr$pos_agnostic <- as_bool(instr$pos_agnostic, default = FALSE)
  instr$keep_all_configurations <- as_bool(
    instr$keep_all_configurations,
    default = FALSE
  )
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

#####

collect_all_layer_configurations <- function(
  GI_data,
  .to_use = c("tech_rep", "bio_rep", "guide_pair"),
  pos_agnostic = FALSE,
  symmetric_analysis_method = "preaverage",
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
      pos_agnostic = pos_agnostic,
      symmetric_analysis_method = symmetric_analysis_method,
      verbose = verbose
    )
  }

  if (length(.collapsable_layers) >= 2) {
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
            message(.n)
          }
          output[[.n]] <- GIScores(
            GI_data,
            collapse_layers = .to_collapse,
            block_layer = .use,
            pos_agnostic = pos_agnostic,
            symmetric_analysis_method = symmetric_analysis_method,
            verbose = verbose
          )
        }
      }
    }
  }
  return(output)
}

#####

find_optimal_configuration <- function(GI_list, keep_all = FALSE) {
  stopifnot(
    "GI_list must contain at least one screen object." = length(GI_list) > 0,
    "GI_list must be named." = !is.null(names(GI_list)) &&
      all(!is.na(names(GI_list))) &&
      all(nzchar(names(GI_list))),
    "GI_list must have unique names." = !anyDuplicated(names(GI_list)),
    "keep_all must be TRUE or FALSE." = is.logical(keep_all) &&
      length(keep_all) == 1 &&
      !is.na(keep_all)
  )

  dupcor <- map_dbl(
    GI_list,
    ~ {
      mean(.x@dupCorrelation, na.rm = TRUE)
    }
  )

  dupcor[is.nan(dupcor)] <- NA_real_

  dupcor_data <- data.table(
    config = names(GI_list),
    dcor = round(dupcor, 3)
  )

  usable <- dupcor[is.finite(dupcor)]

  if (length(usable) == 0) {
    stop(
      "No usable duplicate-correlation estimates found for configuration selection.",
      call. = FALSE
    )
  }

  if (all(usable < 0)) {
    selected_name <- names(usable)[which.max(usable)]
  } else {
    candidates <- usable

    if (any(candidates >= 0)) {
      # if greater or = 0 found, remove <= 0
      candidates <- candidates[candidates >= 0]
    }

    if (any(candidates > 0)) {
      # if any greater 0 found, remove = 0
      candidates <- candidates[candidates > 0]
    }
    # if more than one option remains, choose weakest correlation
    selected_name <- names(candidates)[which.min(candidates)]
  }

  dupcor_data[,
    kept := fifelse(
      config == selected_name,
      "selected",
      ""
    )
  ]

  for (.n in names(GI_list)) {
    GI_list[[.n]]@metadata$dupcor_data <- dupcor_data
  }

  if (isTRUE(keep_all)) {
    return(GI_list)
  } else {
    return(
      GI_list[selected_name]
    )
  }
}

#####
