#' @export

define_screen_type <- function(.checks) {
  
  .type <- data.table::fcase(.checks$gene_sets_equal & 
                               .checks$query_sufficient & 
                               .checks$library_sufficient & 
                               #.attr$checks$stable_library_size & 
                               .checks$sufficient_tests_per_query, "multiplex.symmetric", 
                             !.checks$gene_sets_equal & 
                               .checks$library_sufficient & 
                               .checks$query_sufficient & 
                               #.attr$checks$stable_library_size & 
                               .checks$sufficient_tests_per_query, "multiplex.asymmetric", 
                             !.checks$library_sufficient | 
                               !.checks$stable_library_size | 
                               !.checks$sufficient_tests_per_query# | #avg_tests_per_query <= 50
                             , "fixed", 
                             # ... 
                             default = "unknown")
  
  
  if (.type == "unknown") {
    warning("Unknown screen design! Forcing fixed pair run.")
    .type <- "fixed"
  }
  return(.type)
}

