#' @export

get_screen_attributes <- function(input) {
  
  .a <- list(contrasts = NULL)
  
  #pair <- NULL # to prevent package environment errors
  
  .a$query_genes <- input[, unique(get("query_gene"))]
  .a$library_genes <- input[, unique(get("library_gene"))]
  .a$all_genes <- union(.a$query_genes, .a$library_genes)
  .a$query_genes_not_in_lib <- setdiff(.a$query_genes, .a$library_genes)
  .a$library_genes_not_in_query <- setdiff(.a$library_genes, .a$query_genes)
  
  .a$n_query_genes <- length(.a$query_genes)
  .a$n_lib_genes <- length(.a$library_genes)
  .a$n_all_genes <- length(.a$all_genes)
  
  .a$observations_per_query <- purrr::map_int(purrr::set_names(.a$query_genes), \(.g) {input[query_gene == .g, .N]})
  
  
  .a$all_pairs <- input[, unique(get("pair"))]
  .a$unique_pairs <- input[, unique(sort_gene_pairs(get("query_gene"), get("library_gene")))]
  
  return(.a)
  
}