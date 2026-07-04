fill_gRNA_GIs <- function(gRNA_GI, input, value_var = "GI") {
  #  stopifnot("Value variable not in input data" = value_var %in% colnames(.x),
  #            "Not all dimensions given as structure found in input data!" = all(structure %in% colnames(.x)))

  dims <- c(gRNA_GI@space, gRNA_GI@replicates)

  for (.dim in dims) {
    input[[.dim]] <- factor(input[[.dim]], levels = unique(input[[.dim]]))
  }

  gRNA_GI@data <- acast(
    data = input,
    formula = as.formula(paste0(
      c(gRNA_GI@space, gRNA_GI@replicates),
      collapse = " ~ "
    )),
    value.var = value_var,
    drop = FALSE
  )
  return(gRNA_GI)
}
