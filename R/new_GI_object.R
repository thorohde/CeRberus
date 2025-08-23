#' @export new_GI_object


new_GI_object <- \(input, 
                    contrasts_col = "contrast", 
                    query_col = "query_gene", 
                    lib_col = "library_gene", 
                    bio_rep_col = "bio_rep", 
                    tech_rep_col = "tech_rep", 
                    guide_col = "guide_pair", 
                    gi_col = "GI", 
                    mean_positional_effects = F, 
                    minimal_query_size = 50, 
                    minimal_library_size = 50) {
  
  
  input <- data.table::copy(input)
  
  
  stopifnot("The input object needs to be a data frame." = data.table::is.data.table(input), 
            "The query gene column is not in the provided dataset." = query_col %in% colnames(input), 
            "The library gene column is not in the provided dataset." = lib_col %in% colnames(input))
  
  setnames(input, 
           old = c(contrasts_col, bio_rep_col, tech_rep_col, guide_col, query_col, lib_col, gi_col), 
           new = c("contrast", "bio_rep", "tech_rep", "guide_pair", "query_gene", "library_gene", "GI"), 
           skip_absent = T)
  
  
  .attr <- get_screen_attributes(input)
  
  iwalk(.attr, ~ setattr(input, .y, .x))
  
  input[, replicate := do.call(paste, c(.SD, sep = "_")), 
        .SDcols = names(.attr$rep_layers)]
  
  
  return(input)
}

