#####

compute_dupcor_plot <- function(GI_list, .fpath = NULL, verbose = FALSE) {
  stopifnot(
    "GI_list must contain at least one screen object." = length(GI_list) > 0,
    "verbose must be TRUE or FALSE." = is.logical(verbose) &&
      length(verbose) == 1 &&
      !is.na(verbose)
  )

  for (.n in names(GI_list)) {
    GI_list[[.n]]@metadata$dupcor_plot <- ggplot(
      data = GI_list[[.n]]@metadata$dupcor_data,
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
      plot(GI_list[[.n]]@metadata$dupcor_plot)
    }

    if (!is.null(.fpath)) {
      dir.create(dirname(.fpath), showWarnings = FALSE, recursive = TRUE)

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

#####
