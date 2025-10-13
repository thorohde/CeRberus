#' @export

z_transform <- function(.x, .mean, .sd) {
  
  if (missing(.mean)) {.mean <- mean(.x, na.rm = T)}
  if (missing(.sd)) {.sd <- stats::sd(.x, na.rm = T)}
  
  return((.x - .mean) / .sd)
}
