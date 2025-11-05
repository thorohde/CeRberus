#' @export

find_optimal_configuration <- function(GI_list, verbose = T) {
  
  .n_options_start <- length(GI_list)
  
  .remaining <- GI_list |> map(dupCorrelation) |> map_dbl(mean, na.rm = T)
  
  if (verbose) {
    print(str_c("Start: ", str_c(.remaining, collapse = ", ")))
    }
  
  #if options greater 0 found, remove <= 0
  
  if (any(.remaining >= 0)) {
    .remaining <- .remaining |> keep(.remaining >= 0)
    if (any(.remaining > 0)) {
      .remaining <- .remaining |> keep(.remaining > 0)
    }}
  
  # if more than one option remaining, determine if only one of them has most values above 0
  if (length(.remaining) >= 1) {
    .remaining <- set_names(names(.remaining)) |> 
      map(~ GI_list[[.x]]) |> 
      map(dupCorrelation) |> 
      map(~ {.x[.x >= quantile(.x, 0.1)]}) |>
      map_lgl(~ all(.x > 0))
    
    if (sum(.remaining) == 1) {.remaining <- .remaining[.remaining]}
  } else {
    .remaining <- set_names(names(.remaining)) |> 
      map(~ GI_list[[.x]]) |> 
      map(dupCorrelation) |> 
      map(sum, na.rm = T)
    .remaining <- .remaining[which.max(.remaining)]
  }
  
  if (verbose) {
    print(str_c("End: ", str_c(names(.remaining), collapse = ", ")))
    }
  
  return(GI_list[names(.remaining)])
}