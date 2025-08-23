#' @export flatten_array

flatten_array <- \(x, dnames, value.name = "value") {

  stopifnot("Please provide an array!" = !is.null(dim(x)))

  if (is.null(dimnames(x))) {
    dimnames(x) <- lapply(dim(x), seq_len)
  }

  output <- expand.grid(dimnames(x), KEEP.OUT.ATTRS = F, stringsAsFactors = F)

  output <- data.table::data.table(output)

  output[[value.name]] <- as.vector(x)

  output <- data.table::data.table(data.frame(output))

  if (!missing(dnames)) {
    given_dnames <- length(dnames)
    required_dnames <- length(dim(x))

    if (given_dnames != required_dnames) {
      warning(paste(given_dnames, "dimnames given, but", required_dnames, "dimensions found!"))
    }

    data.table::setnames(output,
                         old = paste0("Var", 1:length(dim(x)))[1:given_dnames],
                         new = dnames)}

return(output)
}
