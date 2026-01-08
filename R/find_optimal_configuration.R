#' @export

find_optimal_configuration <- function(GI_list, verbose = T) {
  
  .x <- GI_list |> map(dupCorrelation) |> map_dbl(mean, na.rm = T)
  
  if (verbose) {
    print("Start: ")
    print(round(.x, 3))
  }
  
  if (any(.x >= 0)) {.x <- keep(.x, .x >= 0)} # ifgreater or = 0 found, remove <= 0
  if (any(.x > 0)) {.x <- .x |> keep(.x > 0)} # if any greater 0 found, remove = 0
  if (length(.x) >= 1) {.x <- .x[which.min(.x)]} # if more than one option remains, choose weakest correlation
  
  if (verbose) {
    print("End: ")
    print(.x)
  }
  
  return(GI_list[names(.x)])
}