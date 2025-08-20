#' @import data.table
#' @importFrom purrr keep map map_dbl map_int set_names reduce
#' @importFrom reshape2 acast
#' @importFrom stats as.formula
#' @importFrom stringr str_c
#' @export make_GI_object


make_GI_object <- \(input, 
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
  


  
  .GI_vals <- if (.run_queried) {c("query_gene", "library_gene", "replicate")} else {c("pair", "replicate")}
  
  if (!is.null(.contrasts)) {.GI_vals <- c(.GI_vals, c("contrast"))}
  
  .rep_layers <- which(.GI_vals == "replicate")


  .GI_vals <- input |> 
    reshape2::acast(formula = as.formula(str_c(.GI_vals, collapse = " ~ ")), 
                    value.var = "GI")
  
  .rep_layers <- replicate_layers(dimnames(.GI_vals)[[.rep_layers]])
  
  ####
  
  if (F) {
    if (!is.null(.contrasts)) {
      .GI_vals <- purrr::map(.contrasts, \(.c) {
        .x <- .GI_vals[,,.c]
        colnames(.x) <- str_c(colnames(.x), "_", .c); .x}) |> 
        purrr::reduce(cbind)
      
      .rep_layers <- purrr::map(.contrasts, \(.c) {
        .x <- .rep_layers
        rownames(.x) <- str_c(rownames(.x), "_", .c)
        .x}) |> purrr::reduce(rbind)
      
    }
  }
  ###
  
  
  
  return(
    list(#attributes = list(query_genes = .query_genes, 
        #                        library_genes = .lib_genes, 
        #                        query_genes_not_in_lib = .query_genes_not_in_lib, 
        #                        library_genes_not_in_query = .lib_genes_not_in_query, 
        #                        all_genes = .all_genes, 
        #                        n_query_genes = .n_query_genes, 
        #                        n_lib_genes = .n_lib_genes, 
        #                        n_all_genes = .n_all_genes, 
        #                        observations_per_query = .observations_per_query, 
        #                        all_pairs = .all_pairs, 
        #                        unique_pairs = .unique_pairs
  #), 
  #replicates = .replicates, 
  replicate_layers = .rep_layers, 
  #contrasts = .contrasts, 
  #checks = .checks, 
  #screen_type = .mode, 
  #run_queried = .run_queried, 
  guide_GIs = .GI_vals
  ))
}

