#' @export new_limma_object

new_limma_object <- function(GI_object) {
  
  .attr <- attributes(GI_object)
  
  if (.attr$screen_type == "unknown") {warning("Unknown screen design! Forcing fixed pair run.")}
  
  output <- c(if (.attr$queried) c("query_gene", "library_gene") else c("pair"), 
              "replicate", 
              if (is.null(.attr$contrasts)) c() else c("contrast")
  )
  
  
  .dim_desc <- map_int(set_names(output), ~ which(output == .x))
  
  output <- GI_object |> 
    reshape2::acast(formula = as.formula(paste0(output, collapse = " ~ ")), 
                    value.var = "GI")
  
  
  data.table::setattr(output, "dim_description", .dim_desc)
  
  walk(c("queried", "contrasts"), ~ data.table::setattr(output, .x, .attr[[.x]]))
  
  data.table::setattr(output, "replicate_layers", replicate_layers(dimnames(output)[[attr(output, "dim_description")[["replicate"]]]]))
  
  data.table::setattr(output, "collapsed_layers", NULL)
  data.table::setattr(output, "block_layer", NULL)

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