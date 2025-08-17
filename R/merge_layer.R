#' @importFrom purrr map set_names
#' @export merge_layer

merge_layer <- function(GI) {
  
  GI <- copy(GI)
  # error handling
  
  .ids <- purrr::map(unique(GI$replicate_layers[,"tech_rep"]), ~ names(which(GI$replicate_layers[,"tech_rep"] == .x)))
  .ids <- .ids |> set_names(str_c("t", 1:length(.ids)))
  
  .new <- empty_array(list(rownames(GI$guide_GIs), names(.ids), 
                           dimnames(GI$guide_GIs)[[3]]))
  
  for (.c in dimnames(.new)[[3]]) {
    for (.tr in names(.ids)) {
      .new[,.tr,.c] <- apply(GI$guide_GIs[,.ids[[.tr]],.c], 1, mean, na.rm = T)}}
  
  GI$guide_GIs <- .new
  
  GI$replicates <- colnames(.new)
  
  GI$replicate_layers <- replicate_layers(GI$replicates)
  
  return(GI)
}