

#' @export get_screen_attributes.R

get_screen_attributes <- function(input, 
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

  stopifnot("The input object needs to be a data frame." = is.data.table(input), 
            "The query gene column is not in the provided dataset." = query_col %in% colnames(input), 
            "The library gene column is not in the provided dataset." = lib_col %in% colnames(input))
  
  setnames(input, 
           old = c(contrasts_col, bio_rep_col, tech_rep_col, guide_col, query_col, lib_col, gi_col), 
           new = c("contrast", "bio_rep", "tech_rep", "guide_pair", "query_gene", "library_gene", "GI"), 
           skip_absent = T)

  .contrasts <- NULL
  #pair <- NULL # to prevent package environment errors
  
  .rep_layers <- set_names(c("contrast", "bio_rep", "tech_rep", "guide_pair")) |>
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
  
  
  ##### 
  
  if (.mode == "unknown") {warning("Unknown screen design! Forcing fixed pair run.")}
  
  .run_queried <- .mode %in% c("symmetric", "asymmetric") & 
    .checks$library_sufficient & 
    .checks$query_sufficient
  
  input[, replicate := apply(.SD, 1, str_c, collapse = "_"), .SDcols = names(.rep_layers)]
  
  
  attr(input, "checks") <- .checks
  
  attr(input, "screen_design") <- list(query_genes = .query_genes, 
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
  
  attr(input, "replicate_design") <- list(#replicates = .replicates, 
    replicate_layers = .rep_layers, 
  contrasts = .contrasts, 

  screen_type = .mode, 
  run_queried = .run_queried#, 
  #guide_GIs = .GI_vals
  )
  
  return(input)
  

  
  
  
}