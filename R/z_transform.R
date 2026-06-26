#' Compute z-scores
#'
#' @description
#' Standardizes a numeric vector by subtracting a mean and dividing by a
#' standard deviation. If `.mean` or `.sd` are not supplied, they are estimated
#' from `.x` with missing values ignored.
#'
#' If the standard deviation is missing or zero, the function returns a vector
#' of `NA_real_` values with the same length as `.x`.
#'
#' @param .x Numeric vector to standardize.
#' @param .mean Optional numeric scalar mean used for centering. If omitted,
#'   `mean(.x, na.rm = TRUE)` is used.
#' @param .sd Optional numeric scalar standard deviation used for scaling. If
#'   omitted, `stats::sd(.x, na.rm = TRUE)` is used.
#'
#' @return A numeric vector of z-scores with the same length as `.x`.
#'
#' @examples
#' z_transform(c(1, 2, 3))
#' z_transform(c(1, 2, 3), .mean = 2, .sd = 1)
#'
#' @export

z_transform <- function(.x, .mean, .sd) {
  if (missing(.mean)) {
    .mean <- mean(.x, na.rm = TRUE)
  }
  if (missing(.sd)) {
    .sd <- stats::sd(.x, na.rm = TRUE)
  }

  if (is.na(.sd) || .sd == 0) {
    return(rep(NA_real_, length(.x)))
  }

  return((.x - .mean) / .sd)
}
