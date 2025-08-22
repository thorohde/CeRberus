#' @importFrom limma duplicateCorrelation
#' @importFrom purrr map
#' 
#' @export compute_dupcor_values

compute_dupcor_values <- function(GI, limit_query = NULL) {
  
  stopifnot("replicate_layers" %in% names(attributes(GI)), 
            is.array(GI), 
            is.array(attr(GI, "replicate_layers")))
  
  
  # case 1: no contrasts, queried
  
  if (is.null(attr(GI, "contrasts")) & attr(GI, "queried")) {
    #print("case1")
    output <- purrr::map(rownames(GI), ~ GI[.x,,])}
  
  # case 2: no contrasts, not queried
  if (is.null(attr(GI, "contrasts")) & !attr(GI, "queried")) {
    #print("case2")
    output <- list(GI)}
  
  # case 3: contrasts, queried
  if (!is.null(attr(GI, "contrasts")) & attr(GI, "queried")) {
    #print("case3")
    output <- purrr::map(rownames(GI), \(.g) {
      purrr::map(dimnames(GI)[[3]], ~ GI[.g,.x])})}
  
  # case 4: contrasts, not queried
  if (!is.null(attr(GI, "contrasts")) & !attr(GI, "queried")) {
    #print("case4")
    output <- purrr::map(GI$contrasts, ~ GI[,,.x])}
  
  
  
  if (!is.null(limit_query) & attr(GI, "queried")) {
    output <- sample(output, min(c(limit_query, length(output))))}
  
  
  #      } else {
  #        purrr::map(sample(rownames(GI), min(c(limit_query, nrow(GI)))), ~ GI[.x,,])
  #      }
  
  
  #}
  
  #print(str(dc))
  
  block_options <- purrr::map(
    purrr::set_names(colnames(attr(GI, "replicate_layers"))), ~ attr(GI, "replicate_layers")[,.x])
  
  
  output <- purrr::map(block_options, \(.b) {
    purrr::map(output, limma::duplicateCorrelation, block = .b, .progress = T) |>
      purrr::map_dbl("consensus.correlation")})
  
  
  
  return(output)
}
