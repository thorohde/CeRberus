#' @export

define_screen_type <- \(.checks) {
  
  data.table::fcase(.checks$gene_sets_equal & 
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
}
