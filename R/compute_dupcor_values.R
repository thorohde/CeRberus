#' @importFrom limma duplicateCorrelation
#' @export compute_dupcor_values

compute_dupcor_values <- function(.data, block_options) {
  
  if (is.array(.data)) {
    .data <- purrr::map(base::rownames(.data), ~ .data[.x,,])
  }
  
  block_options <- purrr::map(
    purrr::set_names(base::colnames(block_options)), ~ block_options[,.x])
  
  output <- purrr::map(block_options, \(.b) {
    purrr::map(.data, limma::duplicateCorrelation, block = .b, .progress = T) |>
      purrr::map_dbl("consensus.correlation")})
  
  base::return(output)}
