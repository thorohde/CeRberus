
#' @importFrom stringr str_c
#' @importFrom purrr map map_dbl
#' @importFrom data.table CJ setnames .SD
#' @export make_GI_object


make_GI_object <- \(.x, 
                    contrasts_col = "contrast", 
                    query_col = "query_gene", 
                    lib_col = "library_gene", 
                    mean_positional_effects = F, 
                    minimal_query_size = 50, 
                    minimal_library_size = 50) {
  
  .x <- copy(.x)
  
  .contrasts <- NULL
  
  stopifnot(query_col %in% colnames(.x), 
            lib_col %in% colnames(.x))
  
  if (contrasts_col %in% colnames(.x)) {
    data.table::setnames(.x, old = contrasts_col, new = "contrast")
    
    .contrasts <- .x[, unique(contrast)]
  }
  
  
  
  data.table::setnames(.x, old = c(query_col, lib_col), new = c("query", "library"))
  
  
  
  
  
  .rep_names <- setdiff(colnames(.x), c("contrast", "query", "library"))
  
  .query_genes <- .x[, base::unique(get("query"))]
  .lib_genes <- .x[, base::unique(get("library"))]
  .all_genes <- base::union(.query_genes, .lib_genes)
  .query_genes_not_in_lib <- setdiff(.query_genes, .lib_genes)
  .lib_genes_not_in_query <- setdiff(.lib_genes, .query_genes)
  
  
  .n_query_genes <- base::length(.query_genes)
  .n_lib_genes <- base::length(.lib_genes)
  .n_all_genes <- base::length(.all_genes)
  
  .observations_per_query <- map_int(set_names(.query_genes), \(.g) {.x[query == .g, .N]})
  
  
  .all_pairs <- .x[, unique(str_c(query, library, sep = ";"))]
  .unique_pairs <- .x[, unique(sort_gene_pairs(query, library))]
  
  .checks <- list(
    gene_sets_equal = base::length(c(.query_genes_not_in_lib, .lib_genes_not_in_query)) == 0,  
    query_sufficient = .n_query_genes >= minimal_query_size, 
    library_sufficient = .n_lib_genes >= minimal_library_size, 
    stable_library_size = sum(.observations_per_query != median(.observations_per_query, na.rm = T)) <= 10, 
    avg_tests_per_query = median(.observations_per_query, .na.rm = T)
    
    #   base::length(.all_pairs) <= 0.9*.n_query_genes*.n_lib_genes
  )
  
  .mode <- fcase(.checks$gene_sets_equal, "symmetric", 
                 !.checks$gene_sets_equal, "asymmetric", 
                 #, "fixed_pairs", 
                 # ... 
                 default = "unknown")
  

  
  if (.mode == "unknown") {warning("Unknown screen design! Forcing fixed pair run.")}
  
  .run_queried <- .mode %in% c("symmetric", "asymmetric") & 
    .checks$library_sufficient & 
    .checks$query_sufficient
  
  
  if ()
  .GI_vals <- list(.query_genes, .lib_genes, .rep_names)
  
  if (!is.null(.contrasts)) {
    .GI_vals <- c(.GI_vals, list(.contrasts))
  }
  
  
  .GI_vals <- base::array(data = NA,
                          dim = purrr::map_dbl(.GI_vals, base::length),
                          dimnames = .GI_vals)
  
  
  .template <- data.table::CJ(query = .query_genes, library = .lib_genes)
  
  if (is.null(.contrasts)) {
    
    .d <- merge(.template, .x, by = c("query", "library"), all.x = T)
    
    for (.r in .rep_names) {
      .GI_vals[,,.r] <- .d[, get(.r)]
      
      if (mean_positional_effects) {
        .GI_vals[,,.r] <- makeSymmetric(.GI_vals[,,.r])
      }
    }
    
  } else {
    
    .x <- purrr::map(set_names(.contrasts), \(.c) {.x[contrast == .c]})
    
    for (.c in .contrasts) {
      .d <- merge(.template, .x[[.c]], by = c("query", "library"), all.x = T)
      
      for (.r in .rep_names) {
        .GI_vals[,,.r,.c] <- .d[, get(.r)]
        
        if (mean_positional_effects) {
          .GI_vals[,,.r,.c] <- makeSymmetric(.GI_vals[,,.r,.c])
        }
      }
    }
  }
  
  
  
  
  return(list(attributes = list(contrasts = .contrasts, 
                                query_genes = .query_genes, 
                                library_genes = .lib_genes, 
                                query_genes_not_in_lib = .query_genes_not_in_lib, 
                                library_genes_not_in_query = .lib_genes_not_in_query, 
                                all_genes = .all_genes, 
                                n_query_genes = .n_query_genes, 
                                n_lib_genes = .n_lib_genes, 
                                n_all_genes = .n_all_genes, 
                                observations_per_query = .observations_per_query, 
                                all_pairs = .all_pairs, 
                                unique_pairs = .unique_pairs, 
                                rep_names = .rep_names
  ), 
  checks = .checks, 
  screen_type = .mode, 
  run_queried = .run_queried, 
  guide_GIs = .GI_vals))
}

