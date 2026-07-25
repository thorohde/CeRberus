#####

compute_dupcor_plot <- function(gi_list, .fpath = NULL, verbose = FALSE) {
  stopifnot(
    "gi_list must contain at least one screen object." = length(gi_list) > 0,
    "verbose must be TRUE or FALSE." = is.logical(verbose) &&
      length(verbose) == 1 &&
      !is.na(verbose)
  )

  for (.n in names(gi_list)) {
    gi_list[[.n]]@metadata$dupcor_plot <- ggplot(
      data = gi_list[[.n]]@metadata$dupcor_data,
      mapping = aes(dcor, config)
    ) +
      theme_light() +
      geom_col(aes(fill = kept)) +
      scale_fill_manual(
        values = purrr::set_names(c("seagreen", "grey80"), c("selected", ""))
      ) +
      geom_vline(xintercept = c(0, 0.25), linetype = "dashed", linewidth = 1) +
      labs(
        #caption = "It is recommended to choose a configuration with most values between 0 and 0.25.",
        x = "Duplicate correlation",
        y = "Limma configuration"
      )

    if (isTRUE(verbose)) {
      plot(gi_list[[.n]]@metadata$dupcor_plot)
    }

    if (!is.null(.fpath)) {
      dir.create(dirname(.fpath), showWarnings = FALSE, recursive = TRUE)

      ggplot2::ggsave(
        filename = .fpath,
        plot = gi_list[[.n]]@metadata$dupcor_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
    }
  }
  return(gi_list)
}

#####

pval_qq_plot <- function(
  gi_obj,
  .fpath = NULL,
  verbose = FALSE,
  ntc = "NTC",
  ntc_string = "target-NTC",
  target_string = "target-target",
  ntc_color = "red",
  target_color = "black",
  title = "QQ-Plot",
  caption = NULL
) {
  stopifnot(
    "gi_obj must be a PosAgnMultiplexScreen object." = methods::is(
      gi_obj,
      "PosAgnMultiplexScreen"
    ),
    "verbose must be TRUE or FALSE." = is.logical(verbose) &&
      length(verbose) == 1 &&
      !is.na(verbose),
    "The symmGeneGIs table must contain a gene_pair column." = "gene_pair" %in%
      colnames(gi_obj@symmGeneGIs),
    "The symmGeneGIs table must contain a pval column." = "pval" %in%
      colnames(gi_obj@symmGeneGIs)
  )

  plot_data <- data.table::copy(gi_obj@symmGeneGIs)

  plot_data <- plot_data[
    !is.na(pval),
    .(
      gene_pair,
      query_gene = stringr::str_split_i(gene_pair, ";", 1),
      library_gene = stringr::str_split_i(gene_pair, ";", 2),
      pval = pmax(as.numeric(pval), .Machine$double.xmin),
      ctrl = data.table::fcase(
        stringr::str_split_i(gene_pair, ";", 1) == ntc | stringr::str_split_i(gene_pair, ";", 2) == ntc ,
        ntc_string                                                                                      ,
        default = target_string
      )
    )
  ]

  stopifnot(
    "symmGeneGIs must contain at least one non-missing p-value." = nrow(
      plot_data
    ) >
      0
  )

  plot_data <- plot_data[order(ctrl, pval, decreasing = FALSE)]
  plot_data <- plot_data[,
    .(
      gene_pair,
      pval,
      obs = pval,
      exp = stats::ppoints(.N)
    ),
    by = .(ctrl)
  ]

  inflation_summary <- plot_data[,
    .(
      n = .N,
      lambda = stats::qchisq(stats::median(1 - pval), df = 1) /
        stats::qchisq(0.5, df = 1),
      min_p = min(pval),
      median_p = stats::median(pval)
    ),
    by = .(ctrl)
  ]

  inflation_caption <- paste0(
    inflation_summary$ctrl,
    ": n=",
    inflation_summary$n,
    ", lambda=",
    formatC(inflation_summary$lambda, digits = 3, format = "f"),
    ", median p=",
    formatC(inflation_summary$median_p, digits = 3, format = "f")
  )
  inflation_caption <- paste(inflation_caption, collapse = " | ")

  full_caption <- if (is.null(caption) || identical(caption, "")) {
    inflation_caption
  } else {
    paste(caption, inflation_caption, sep = "\n")
  }

  dot_colors <- purrr::set_names(
    c(ntc_color, target_color),
    c(ntc_string, target_string)
  )

  gi_obj@metadata$qq_plot_data <- plot_data
  gi_obj@metadata$qq_inflation_summary <- inflation_summary

  gi_obj@metadata$qq_plot <- ggplot2::ggplot(
    data = plot_data,
    mapping = ggplot2::aes(
      x = -log10(exp),
      y = -log10(obs),
      color = ctrl
    )
  ) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_color_manual(values = dot_colors) +
    ggplot2::labs(
      title = title,
      caption = full_caption,
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      color = NULL
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue") +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "bottom")

  if (isTRUE(verbose)) {
    plot(gi_obj@metadata$qq_plot)
  }

  if (!is.null(.fpath)) {
    dir.create(dirname(.fpath), showWarnings = FALSE, recursive = TRUE)

    ggplot2::ggsave(
      filename = .fpath,
      plot = gi_obj@metadata$qq_plot,
      width = 8,
      height = 5,
      dpi = 300
    )
  }

  return(gi_obj)
}

#####
