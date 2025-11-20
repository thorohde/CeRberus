#' @export

collect_all_layer_configurations <- function(GI_data, 
                                    .collapsable_layers = c("tech_rep", "bio_rep", "guide_pair"), 
                                    .to_use = c("tech_rep", "bio_rep", "guide_pair"), 
                                    make_pos_agnostic = F) {
  
  .collapsable_layers <- intersect(.collapsable_layers, colnames(GI_data))
  .to_use = intersect(.to_use, colnames(GI_data))
  
  output <- list()
  
  for (.use in .to_use) {
    
    output[[str_c("default_", .use, "_used")]] <- GIScores(
      data.table::copy(GI_data), 
      block_layer = .use, 
      pos_agnostic = make_pos_agnostic)
  }
  
  for (.i in seq_len(length(.collapsable_layers)-1)) {
    for (.to_collapse in combn(.collapsable_layers, .i, simplify = F)) {
      for (.use in setdiff(.collapsable_layers, .to_collapse)) {
        output[[str_c(str_c(.to_collapse, collapse = "_"), "_collapsed_", .use, "_used")]] <- GIScores(
          data.table::copy(GI_data), 
          collapse_layers = .to_collapse, 
          block_layer = .use, 
          pos_agnostic = make_pos_agnostic)
      }}}
  return(output)
  }