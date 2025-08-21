#' @importFrom purrr
#' @importFrom reshape2 acast
#' @importFrom stats as.formula
#' @export build_limma_object

build_limma_object <- function(GI_object) {
  
  .attr <- attributes(GI_object)
  
  if (.attr$screen_type == "unknown") {warning("Unknown screen design! Forcing fixed pair run.")}
  
  output <- if (.attr$queried) {
    c("query_gene", "library_gene", "replicate")
  } else {c("pair", "replicate")}
  
  
  if (!is.null(.attr$contrasts)) {
    output <- c(output, c("contrast"))
    }
  
  .dim_desc <- map_int(set_names(output), ~ which(output == .x))
  
  output <- GI_object |> 
    reshape2::acast(formula = as.formula(paste0(output, collapse = " ~ ")), 
                    value.var = "GI")
  
  setattr(output, "dim_description", .dim_desc)
          
  attr(output, "replicate_layers") <- replicate_layers(dimnames(output)[[attr(output, "dim_description")[["replicate"]]]])
  
  ####
  
  if (F) {
    if (!is.null(.attr$contrasts)) {
      output <- purrr::map(.attr$contrasts, \(.c) {
        .x <- output[,,.c]
        colnames(.x) <- paste0(colnames(.x), "_", .c); .x}) |> 
        purrr::reduce(cbind)
      
      .rep_layers <- purrr::map(.attr$contrasts, \(.c) {
        .x <- .rep_layers
        rownames(.x) <- paste0(rownames(.x), "_", .c)
        .x}) |> purrr::reduce(rbind)
    }
  }
  ###
  
  
  
  
  return(output)
  
}