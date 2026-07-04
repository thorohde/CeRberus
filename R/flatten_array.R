#' Flatten an array to long format
#'
#' @description
#' Converts an array into a long-format `data.table` with one row per array
#' cell. Array dimension names become columns and array values are stored in a
#' separate value column.
#'
#' If the input array has no `dimnames`, integer sequence labels are generated
#' for each dimension before flattening.
#'
#' @param x An array or array-like object with a non-`NULL` `dim` attribute.
#' @param dnames Optional character vector used to rename the generated
#'   dimension columns. Its length should match the number of dimensions in
#'   `x`.
#' @param value.name Character scalar naming the output column that stores array
#'   values. Defaults to `"value"`.
#'
#' @return A `data.table` in long format containing one column per array
#'   dimension and one value column named by `value.name`.
#'
#' @examples
#' x <- array(
#'   1:4,
#'   dim = c(2, 2),
#'   dimnames = list(c("A", "B"), c("C", "D"))
#' )
#' flatten_array(x, dnames = c("query_gene", "library_gene"))
#'
#' @export

flatten_array <- \(x, dnames, value.name = "value") {
  stopifnot("Please provide an array!" = !is.null(dim(x)))

  if (is.null(dimnames(x))) {
    dimnames(x) <- lapply(dim(x), seq_len)
  }

  required_dnames <- length(dim(x))

  if (required_dnames == 2) {
    # Force 2D matrix to melt cleanly into a long data.table
    output <- data.table::as.data.table(as.data.frame.table(
      x,
      responseName = value.name
    ))
  } else {
    # 3D+ arrays use the native data.table array melting
    output <- data.table::as.data.table(x, value.name = value.name)
  }

  if (!missing(dnames)) {
    given_dnames <- length(dnames)

    if (given_dnames != required_dnames) {
      warning(paste(
        given_dnames,
        "dimnames given, but",
        required_dnames,
        "dimensions found!"
      ))
    }

    rename_n <- min(given_dnames, required_dnames)
    dim_cols <- names(output)[seq_len(rename_n)]

    data.table::setnames(
      output,
      old = dim_cols,
      new = dnames[seq_len(rename_n)]
    )
  }

  return(output)
}
