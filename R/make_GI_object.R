#' @import data.table
#' @importFrom purrr map map_dbl map_int set_names reduce
#' @importFrom reshape2 acast
#' @importFrom stats as.formula
#' @importFrom stringr str_c
#' @export make_GI_object


make_GI_object <- \(.x, 
                    contrasts_col = "contrast", 
                    query_col = "query_gene", 
                    lib_col = "library_gene", 
                    mean_positional_effects = F, 
                    minimal_query_size = 50, 
                    minimal_library_size = 50) {
  
  .x <- data.table::copy(.x)
  
  .contrasts <- NULL
  #pair <- NULL # to prevent package environment errors
  
  
  stopifnot(query_col %in% colnames(.x), 
            lib_col %in% colnames(.x))
  
  if (contrasts_col %in% colnames(.x)) {
    data.table::setnames(.x, old = contrasts_col, new = "contrast")
    
    .contrasts <- .x[, unique(get("contrast"))]
  }
  
  
  
  data.table::setnames(.x, old = c(query_col, lib_col), new = c("query", "library"))
  .x[, pair := stringr::str_c(get("query"), ";", get("library"))]
  
  .replicates <- setdiff(colnames(.x), c("contrast", "query", "library", "pair"))
  
  .query_genes <- .x[, unique(get("query"))]
  .lib_genes <- .x[, unique(get("library"))]
  .all_genes <- union(.query_genes, .lib_genes)
  .query_genes_not_in_lib <- setdiff(.query_genes, .lib_genes)
  .lib_genes_not_in_query <- setdiff(.lib_genes, .query_genes)
  
  
  .n_query_genes <- length(.query_genes)
  .n_lib_genes <- length(.lib_genes)
  .n_all_genes <- length(.all_genes)
  
  .observations_per_query <- purrr::map_int(purrr::set_names(.query_genes), \(.g) {.x[query == .g, .N]})
  
  
  .all_pairs <- .x[, unique(get("pair"))]
  .unique_pairs <- .x[, unique(sort_gene_pairs(get("query"), get("library")))]
  
  .checks <- list(
    gene_sets_equal = length(c(.query_genes_not_in_lib, .lib_genes_not_in_query)) == 0,  
    query_sufficient = .n_query_genes >= minimal_query_size, 
    library_sufficient = .n_lib_genes >= minimal_library_size, 
    stable_library_size = sum(.observations_per_query != stats::median(.observations_per_query, na.rm = T)) <= 10, 
    sufficient_tests_per_query = sum(.observations_per_query >= minimal_library_size) >= 0.95 * length(.observations_per_query), 
    
    
    avg_tests_per_query = stats::median(.observations_per_query, .na.rm = T)
    
    
    #   length(.all_pairs) <= 0.9*.n_query_genes*.n_lib_genes
  )
  
  .mode <- data.table::fcase(.checks$gene_sets_equal & 
                               .checks$query_sufficient & 
                               .checks$library_sufficient & 
                               .checks$stable_library_size & 
                               .checks$sufficient_tests_per_query, "symmetric", 
                             !.checks$gene_sets_equal & 
                               .checks$library_sufficient & 
                               .checks$query_sufficient & 
                               .checks$stable_library_size & 
                               .checks$sufficient_tests_per_query, "asymmetric", 
                             !.checks$library_sufficient | 
                               !.checks$stable_library_size | 
                               !.checks$sufficient_tests_per_query# | #avg_tests_per_query <= 50
                             , "fixed_pairs", 
                             # ... 
                             default = "unknown")
  
  
  
  if (.mode == "unknown") {warning("Unknown screen design! Forcing fixed pair run.")}
  
  .run_queried <- .mode %in% c("symmetric", "asymmetric") & 
    .checks$library_sufficient & 
    .checks$query_sufficient
  
  .GI_vals <- if (.run_queried) {c("query", "library", "replicate")} else {c("pair", "replicate")}
  
  
  if (!is.null(.contrasts)) {.GI_vals <- c(.GI_vals, c("contrast"))}
  
  .GI_vals <- as.formula(str_c(.GI_vals, collapse = " ~ "))
  
  .GI_vals <- .x |> melt.data.table(measure.vars = .replicates, 
                                    variable.name = "replicate", 
                                    value.name = "GI", 
                                    variable.factor = F, 
                                    value.factor = F) |> 
    reshape2::acast(formula = .GI_vals, 
                    value.var = "GI")
  
  
  
  .rep_layers <- replicate_layers(.replicates)
  
  
  
  if (F) {
    if (!is.null(.contrasts)) {
      .GI_vals <- purrr::map(.contrasts, \(.c) {
        .x <- .GI_vals[,,.c]
        colnames(.x) <- str_c(colnames(.x), "_", .c); .x}) %>% 
        purrr::reduce(cbind)
      
      .rep_layers <- purrr::map(.contrasts, \(.c) {
        .x <- .rep_layers
        rownames(.x) <- str_c(rownames(.x), "_", .c)
        .x}) %>% purrr::reduce(rbind)
      
    }
  }
  
  
  
  
  return(list(attributes = list(query_genes = .query_genes, 
                                library_genes = .lib_genes, 
                                query_genes_not_in_lib = .query_genes_not_in_lib, 
                                library_genes_not_in_query = .lib_genes_not_in_query, 
                                all_genes = .all_genes, 
                                n_query_genes = .n_query_genes, 
                                n_lib_genes = .n_lib_genes, 
                                n_all_genes = .n_all_genes, 
                                observations_per_query = .observations_per_query, 
                                all_pairs = .all_pairs, 
                                unique_pairs = .unique_pairs
  ), 
  replicates = .replicates, 
  replicate_layers = .rep_layers, 
  contrasts = .contrasts, 
  checks = .checks, 
  screen_type = .mode, 
  run_queried = .run_queried, 
  guide_GIs = .GI_vals))
}

