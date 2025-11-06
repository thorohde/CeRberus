#' @export

abba_cor <- \(.x) {
  stopifnot(length(dim(.x)) == 2 & dim(.x)[1] == dim(.x)[2])
  return(cor(as.vector(.x), as.vector(t(.x)), use = "pairwise.complete.obs"))
}
