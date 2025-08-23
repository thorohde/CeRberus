#' @export compute_dupcor_values

compute_dupcor_values <- function(GI, sample_query = NULL) {
  
  stopifnot("replicate_layers" %in% names(attributes(GI)), 
            is.array(GI), 
            is.array(attr(GI, "replicate_layers")))
  
  if (is.null(attr(GI, "contrasts"))) {
    
    if (attr(GI, "queried")) {
      output <- purrr::map(rownames(GI), ~ GI[.x,,])
    } else {
      output <- list(GI)
    }
  }
  
  if (!is.null(attr(GI, "contrasts"))) {
    
    if (attr(GI, "queried")) {
      output <- purrr::map(rownames(GI), \(.g) {
        purrr::map(dimnames(GI)[[3]], ~ GI[.g,.x])})
    } else {
      output <- purrr::map(GI$contrasts, ~ GI[,,.x])
    }
  }
  
  if (!is.null(sample_query) & attr(GI, "queried")) {
    output <- sample(output, min(c(sample_query, length(output))))}
  
  
  block_options <- if (is.null(attr(GI, "block_layer"))) {
    colnames(attr(GI, "replicate_layers"))
  } else {
    attr(GI, "block_layer")
  }
  
  block_options <- purrr::map(purrr::set_names(block_options), ~ attr(GI, "replicate_layers")[,.x])
  
  
  
  output <- purrr::map(block_options, \(.b) {
    purrr::map(output, limma::duplicateCorrelation, block = .b) |>
      purrr::map_dbl("consensus.correlation")})
  
  
  
  return(output)
}
