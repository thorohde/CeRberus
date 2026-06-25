#' @export

z_transform <- function(.x, .mean, .sd) {
  if (missing(.mean)) {
    .mean <- mean(.x, na.rm = T)
  }
  if (missing(.sd)) {
    .sd <- stats::sd(.x, na.rm = T)
  }

  if (is.na(.sd) || .sd == 0) {
    return(rep(NA_real_, length(.x)))
  }

  return((.x - .mean) / .sd)
}
