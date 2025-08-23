#' @export collapse_layer


collapse_layer <- function(input, .collapse) {
  
  input <- data.table::copy(input)
  
  stopifnot(.collapse %in% colnames(attr(input, "replicate_layers")))
  
  .old_reps <- colnames(attr(input, "replicate_layers"))
  .new_reps <- setdiff(.old_reps, .collapse)
  
  output <- input |> flatten_array(names(attr(input, "dim_description")), "GI")
  
  output[, c(.old_reps) := data.table::tstrsplit(replicate, "_")]
  
  output[, replicate := do.call(paste, c(data.table::.SD, sep = "_")), .SDcols = .new_reps]
  
  output <- output |> 
    acast(formula = as.formula("query_gene ~ library_gene ~ replicate"), 
          value.var = "GI", 
          fun.aggregate = mean)
  
  walk(c("block_layer", "collapsed_layers", "queried", "dim_description"), ~ data.table::setattr(output, .x, attr(input, .x)))
  
  data.table::setattr(output, "collapsed_layers", c(attr(output, "collapsed_layers"), .collapse))

  data.table::setattr(output, "replicate_layers", replicate_layers(dimnames(output)[[attr(output, "dim_description")[["replicate"]]]]))
  
  return(output)
}