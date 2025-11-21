#' @export

usable_for_limma <- \(.x) {
  
  .usable <- is.array(.x) & sum(is.na(.x)) / length(.x) <= 0.9
  
  return(.usable)
  }
