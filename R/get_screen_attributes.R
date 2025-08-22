

#' @export get_screen_attributes

get_screen_attributes <- function(input, 
                                  #mean_positional_effects = F, 
                                  minimal_query_size = 50, 
                                  minimal_library_size = 50) {
  
  
  
  .attr <- list(contrasts = NULL)
  
  
  
  #pair <- NULL # to prevent package environment errors
  
  .attr$rep_layers <- set_names(c("contrast", "bio_rep", "tech_rep", "guide_pair")) |>
    keep(~ .x %in% colnames(input)) |>
    map(~ input[, sort(unique(get(.x)))])
  
  input[, pair := stringr::str_c(get("query_gene"), ";", get("library_gene"))]
  
  .query_genes <- input[, unique(get("query_gene"))]
  .lib_genes <- input[, unique(get("library_gene"))]
  .all_genes <- union(.query_genes, .lib_genes)
  .query_genes_not_in_lib <- setdiff(.query_genes, .lib_genes)
  .lib_genes_not_in_query <- setdiff(.lib_genes, .query_genes)
  
  
  .n_query_genes <- length(.query_genes)
  .n_lib_genes <- length(.lib_genes)
  .n_all_genes <- length(.all_genes)
  
  .observations_per_query <- purrr::map_int(purrr::set_names(.query_genes), \(.g) {input[query_gene == .g, .N]})
  
  
  .all_pairs <- input[, unique(get("pair"))]
  .unique_pairs <- input[, unique(sort_gene_pairs(get("query_gene"), get("library_gene")))]
  
  .attr$checks <- list(
    gene_sets_equal = (length(.query_genes_not_in_lib) <= 0.02*.n_query_genes) & 
      (length(.lib_genes_not_in_query) <= 0.02*.n_lib_genes), 
    query_sufficient = .n_query_genes >= minimal_query_size, 
    library_sufficient = .n_lib_genes >= minimal_library_size, 
    stable_library_size = sum(.observations_per_query != stats::median(.observations_per_query, na.rm = T)) <= 10, 
    sufficient_tests_per_query = sum(.observations_per_query >= minimal_library_size) >= 0.95 * length(.observations_per_query), 
    
    
    avg_tests_per_query = stats::median(.observations_per_query, .na.rm = T)
    
    
    #   length(.all_pairs) <= 0.9*.n_query_genes*.n_lib_genes
  )
  
  .attr$screen_type <- data.table::fcase(.attr$checks$gene_sets_equal & 
                                           .attr$checks$query_sufficient & 
                                           .attr$checks$library_sufficient & 
                                           #.attr$checks$stable_library_size & 
                                           .attr$checks$sufficient_tests_per_query, "symmetric", 
                                         !.attr$checks$gene_sets_equal & 
                                           .attr$checks$library_sufficient & 
                                           .attr$checks$query_sufficient & 
                                           #.attr$checks$stable_library_size & 
                                           .attr$checks$sufficient_tests_per_query, "asymmetric", 
                                         !.attr$checks$library_sufficient | 
                                           !.attr$checks$stable_library_size | 
                                           !.attr$checks$sufficient_tests_per_query# | #avg_tests_per_query <= 50
                                         , "fixed_pairs", 
                                         # ... 
                                         default = "unknown")
  
  
  ##### 
  
  
  
  .attr$queried <- .attr$screen_type %in% c("symmetric", "asymmetric") & 
    .attr$checks$library_sufficient & 
    .attr$checks$query_sufficient
  
  .attr$screen_design <- list(query_genes = .query_genes, 
                              library_genes = .lib_genes, 
                              query_genes_not_in_lib = .query_genes_not_in_lib, 
                              library_genes_not_in_query = .lib_genes_not_in_query, 
                              all_genes = .all_genes, 
                              n_query_genes = .n_query_genes, 
                              n_lib_genes = .n_lib_genes, 
                              n_all_genes = .n_all_genes, 
                              observations_per_query = .observations_per_query, 
                              all_pairs = .all_pairs, 
                              unique_pairs = .unique_pairs)
  
  return(.attr)
  
}