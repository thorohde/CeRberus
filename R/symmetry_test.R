#' @export

symmetry_test <- \(GI_obj, cutoff = 0.99) {
  
  stopifnot(screenType(GI_obj) == "multiplex.symmetric")
  
  
  .test <- map_lgl(set_names(replicates(GI_obj)), \(.r) {
    .x <- guideGIs(GI_obj)[,,.r]
    if (all(.x == t(.x)) || all(dplyr::near(.x, t(.x)))) {return(T)} else {
      return(cor(as.vector(.x), as.vector(t(.x)), use = "pairwise.complete.obs") >= cutoff)}})

return(all(.test))  

}
