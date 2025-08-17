#' @importFrom abind abind
#' @export makeSymmetric

makeSymmetric <- function(.x) {
  base::apply(
    abind::abind(.x, base::t(.x), along = 3), 1:2, 
    base::mean, na.rm = T)
}

