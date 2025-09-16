#' @export

GIScores <- function(input, 
                     contrasts_col = "contrast", 
                     query_col = "query_gene", 
                     lib_col = "library_gene", 
                     bio_rep_col = "bio_rep", 
                     tech_rep_col = "tech_rep", 
                     guide_col = "guide_pair", 
                     gi_col = "GI", 
                     mean_positional_effects = F, 
                     min_query_size = 50, 
                     min_library_size = 50) {
  
  input <- data.table::copy(input)
  
  input[, pair := paste0(get("query_gene"), ";", get("library_gene"))]  
  
  input[, replicate := do.call(paste, c(.SD, sep = "_")), 
        .SDcols = intersect(c("contrast", "bio_rep", "tech_rep", "guide_pair"), colnames(input))]
  
  
  stopifnot("The input object needs to be a data frame." = data.table::is.data.table(input), 
            "The query gene column is not in the provided dataset." = query_col %in% colnames(input), 
            "The library gene column is not in the provided dataset." = lib_col %in% colnames(input))
  
  
  setnames(input, 
           old = c(contrasts_col, bio_rep_col, tech_rep_col, guide_col, query_col, lib_col, gi_col), 
           new = c("contrast", "bio_rep", "tech_rep", "guide_pair", "query_gene", "library_gene", "GI"), 
           skip_absent = T)
  
  
  .a <- get_screen_attributes(input)
  
  
  
  .checks <- list(
    gene_sets_equal = (length(.a$query_genes_not_in_lib) <= 0.02*.a$n_query_genes) & 
      (length(.a$lib_genes_not_in_query) <= 0.02*.a$n_lib_genes), 
    query_sufficient = .a$n_query_genes >= min_query_size, 
    library_sufficient = .a$n_lib_genes >= min_library_size, 
    stable_library_size = sum(.a$observations_per_query != stats::median(.a$observations_per_query, na.rm = T)) <= 10, 
    sufficient_tests_per_query = sum(.a$observations_per_query >= min_library_size) >= 0.95 * length(.a$observations_per_query), 
    avg_tests_per_query = stats::median(.a$observations_per_query, .na.rm = T)
    
    #   length(.all_pairs) <= 0.9*.n_query_genes*.n_lib_genes
  )
  
  
  .type <- define_screen_type(.checks)
  
  if (.type == "unknown") {
    warning("Unknown screen design! Forcing fixed pair run.")
    .type <- "fixed"
  }
  
  
  GI_object <- new("GIScores")
  
  screen_attributes(GI_object) <- .a
  checks(GI_object) <- .checks
  screen_type(GI_object) <- .type
  
  .f <- c(if (grepl("multiplex", .type)) c("query_gene", "library_gene") else c("pair"), 
          "replicate"#, 
          #if (is.null(screen_attributes(GI_object)$contrasts)) c() else c("contrast")
  )
  
  dim_description(GI_object) <- purrr::map_int(set_names(.f), ~ which(.f == .x))
  
  guideGIs(GI_object) <- input |> 
    reshape2::acast(formula = as.formula(paste0(.f, collapse = " ~ ")), 
                    value.var = "GI")
  
  layers(GI_object) <- replicate_layers(dimnames(guideGIs(GI_object))[[dim_description(GI_object)[["replicate"]]]])
  
  ###
  #data.table::setattr(GI_object, "collapsed_layers", NULL)
  #data.table::setattr(GI_object, "block_layer", NULL)
  ###
  
  return(GI_object)
}
