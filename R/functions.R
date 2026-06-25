#' @export

normalizeReadcounts <- function(readcounts, cf1 = 100, cf2 = 1) {
  # replaced 1e6 with 100 x length(readcounts), cf1 = 1e6, cf2 = 0.5
  x <- log2(
    (readcounts / sum(readcounts, na.rm = TRUE)) *
      cf1 *
      length(readcounts) +
      cf2
  ) #NA will stay NA
  if (!all(is.na(x))) {
    #not run if replicate is missing (= all NA)
    x <- x - min(x, na.rm = TRUE)
  } #smallest value is 0 regardless of cf2
  x
}


remove_PCs <- \(.x, to_remove = NA, .center = TRUE, .scale = TRUE) {
  pca_result <- prcomp(.x, center = .center, scale. = .scale)
  pcs <- pca_result$x
  center <- pca_result$center
  scale <- pca_result$scale

  if (all(is.na(to_remove))) {
    return(.x)
  }

  to_remove <- unique(to_remove) # Validate to_remove indices
  if (any(to_remove > ncol(pcs))) {
    stop(
      "Some PCs to remove exceed the number of available principal components."
    )
  }

  pcs[, to_remove] <- 0
  .x <- pcs %*% t(pca_result$rotation)
  if (.scale) {
    .x <- sweep(.x, 2, scale, FUN = "*")
  }
  if (.center) {
    .x <- sweep(.x, 2, center, FUN = "+")
  }
  return(.x)
}


read_instructions <- function(yaml_fpath) {
  stopifnot("Instruction file does not exist!" = file.exists(yaml_fpath))

  instr <- read_yaml(file = yaml_fpath)

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
    "Given scores file does not exist!" = file.exists(instr$scores_file),
    "Output directory required!" = "output_directory" %in% names(instr),
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

  if (endsWith(instr$output_directory, "/")) {
    warning("Please remove the trailing '/' in the output_directory path!")
    instr$output_directory <- gsub("/$", "", instr$output_directory)
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
      config %in% names(.x) , "*" ,
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


compute_dupcor_plot <- function(GI_list, .fpath = NULL, verbose = FALSE) {
  for (.n in names(GI_list)) {
    GI_list[[.n]]@metadata$dupcor_plot <- ggplot(
      data = GI_list[[.n]]@metadata$dupcor_data,
      mapping = aes(dcor, config)
    ) +
      theme_light() +
      geom_col(aes(fill = kept)) +
      scale_fill_manual(values = c("TRUE" = "seagreen", "FALSE" = "NA")) +
      geom_vline(xintercept = c(0, 0.25), linetype = "dashed", linewidth = 1) +
      labs(
        #caption = "It is recommended to choose a configuration with most values between 0 and 0.25.",
        x = "Duplicate correlation",
        y = "Limma configuration"
      )

    if (isTRUE(verbose)) {
      plot(GI_list[[.n]]@metadata$dupcor_plot)
    }

    if (!is.null(.fpath)) {
      dir.create(dirname(.fpath), F, T)

      ggplot2::ggsave(
        filename = .fpath,
        plot = GI_list[[.n]]@metadata$dupcor_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
    }
  }
  return(GI_list)
}
