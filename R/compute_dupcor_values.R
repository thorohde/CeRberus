#' @importFrom limma duplicateCorrelation
#' @importFrom purrr map
#' 
#' @export compute_dupcor_values

compute_dupcor_values <- function(GI, limit_query = NULL) {
  
  stopifnot("guide_GIs" %in% names(GI), 
            "replicate_layers" %in% names(GI), 
            "contrasts" %in% names(GI), 
            is.array(GI$guide_GIs), 
            is.array(GI$replicate_layers))
  
  
  # case 1: no contrasts, queried
  
  if (is.null(GI$contrasts) & GI$run_queried) {
    #print("case1")
    output <- purrr::map(rownames(GI$guide_GIs), ~ GI$guide_GIs[.x,,])}
  
  # case 2: no contrasts, not queried
  if (is.null(GI$contrasts) & !GI$run_queried) {
    #print("case2")
    output <- list(GI$guide_GIs)}
  
  # case 3: contrasts, queried
  if (!is.null(GI$contrasts) & GI$run_queried) {
    #print("case3")
    output <- purrr::map(rownames(GI$guide_GIs), \(.g) {
      purrr::map(dimnames(GI$guide_GIs)[[3]], ~ GI$guide_GIs[.g,.x])})}
  
  # case 4: contrasts, not queried
  if (!is.null(contrasts) & !GI$run_queried) {
    #print("case4")
    output <- purrr::map(GI$contrasts, ~ GI$guide_GIs[,,.x])}
  
  
  
  if (!is.null(limit_query) & GI$run_queried) {
    output <- sample(output, min(c(limit_query, length(output))))}
  
  
  #      } else {
  #        purrr::map(sample(rownames(GI$guide_GIs), min(c(limit_query, nrow(GI$guide_GIs)))), ~ GI$guide_GIs[.x,,])
  #      }
  
  
  #}
  
  #print(str(dc))
  
  block_options <- purrr::map(
    purrr::set_names(colnames(GI$replicate_layers)), ~ GI$replicate_layers[,.x])
  
  
  output <- purrr::map(block_options, \(.b) {
    purrr::map(output, limma::duplicateCorrelation, block = .b, .progress = T) |>
      purrr::map_dbl("consensus.correlation")})
  
  
  
  return(output)
}
