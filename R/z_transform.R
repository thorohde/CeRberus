#' Compute z-scores
#'
#' @description
#' Standardizes a numeric vector by subtracting a mean and dividing by a
#' standard deviation. If `.mean` or `.sd` are not supplied, they are estimated
#' from `.x` with missing values ignored. When `outlier_quantile` is supplied,
#' the standard deviation is instead estimated from values within the specified
#' central quantile interval.
#'
#' If the standard deviation is missing or zero, the function returns a vector
#' of `NA_real_` values with the same length as `.x`.
#'
#' @param .x Numeric vector to standardize.
#' @param .mean Optional numeric scalar mean used for centering. If omitted,
#'   `mean(.x, na.rm = TRUE)` is used.
#' @param .sd Optional numeric scalar standard deviation used for scaling. If
#'   omitted, `stats::sd(.x, na.rm = TRUE)` is used.
#' @param outlier_quantile Optional numeric scalar in the interval `[0, 0.5)`.
#'   When supplied, the standard deviation is estimated from values between the
#'   `outlier_quantile` and `1 - outlier_quantile` quantiles of `.x`, overriding
#'   `.sd`. Missing values are ignored.
#'
#' @return A numeric vector of z-scores with the same length as `.x`.
#'
#' @examples
#' z_transform(c(1, 2, 3))
#' z_transform(c(1, 2, 3), .mean = 2, .sd = 1)
#'
#' @export

#####

z_transform <- function(.x, .mean, .sd, outlier_quantile = NULL) {
  if (missing(.mean)) {
    .mean <- mean(.x, na.rm = TRUE)
  }
  if (missing(.sd)) {
    .sd <- stats::sd(.x, na.rm = TRUE)
  }

  if (!is.null(outlier_quantile)) {
    if (
      length(outlier_quantile) != 1L ||
        !is.numeric(outlier_quantile) ||
        !is.finite(outlier_quantile) ||
        outlier_quantile < 0 ||
        outlier_quantile >= 0.5
    ) {
      stop("`outlier_quantile` must be one finite numeric value in [0, 0.5).")
    }

    quantile_bounds <- stats::quantile(
      .x,
      probs = c(outlier_quantile, 1 - outlier_quantile),
      na.rm = TRUE
    )
    central_values <- .x[data.table::between(
      .x,
      quantile_bounds[[1L]],
      quantile_bounds[[2L]]
    )]
    .sd <- stats::sd(central_values, na.rm = TRUE)
  }

  if (is.na(.sd) || .sd == 0) {
    return(rep(NA_real_, length(.x)))
  }

  return((.x - .mean) / .sd)
}
