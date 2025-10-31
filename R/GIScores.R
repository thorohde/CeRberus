#' @export

GIScores <- function(input, 
                     query_col = "query_gene", 
                     lib_col = "library_gene", 
                     bio_rep_col = "bio_rep", 
                     tech_rep_col = "tech_rep", 
                     guide_col = "guide_pair", 
                     gi_col = "GI", 
                     collapse_layers = NULL, 
                     block_layer = NULL, 
                     #mean_positional_effects = F, 
                     min_query_size = 50, 
                     min_library_size = 50) {
  
  #message(paste0("Collapsing layers: ", collapse_layers, collapse = ", "))
  #message(paste0(block_layer, " used as blocks"))
  
  input <- data.table::copy(input)
  

  
  stopifnot("The input object needs to be a data frame." = data.table::is.data.table(input), 
            "The query gene column is not in the provided dataset." = query_col %in% colnames(input), 
            "The library gene column is not in the provided dataset." = lib_col %in% colnames(input))
  
  setnames(input, 
           old = c(bio_rep_col, tech_rep_col, guide_col, query_col, lib_col, gi_col), 
           new = c("bio_rep", "tech_rep", "guide_pair", "query_gene", "library_gene", "GI"), 
           skip_absent = T)
  
  
  input[, pair := paste0(get("query_gene"), ";", get("library_gene"))]
  
  #####
  

  if (!is.null(collapse_layers)) {
    input <- input[, .SD, .SDcols = setdiff(colnames(input), collapse_layers)]
    input <- input[, .(GI = mean(GI, na.rm = T)), by = setdiff(colnames(input), "GI")]
  }
  
  

  
  ######
  
  input[, replicate := do.call(paste, c(.SD, sep = "_")), 
        .SDcols = intersect(c("bio_rep", "tech_rep", "guide_pair"), colnames(input))]
  

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
  
  GI_obj <- new("ScreenBase")
  
  screen_attributes(GI_obj) <- .a
  checks(GI_obj) <- .checks
  screenType(GI_obj) <- .type
  
  structure(GI_obj) <- c(if (grepl("multiplex", .type)) c("query_gene", "library_gene") else c("pair"), "replicate")
  
  guideGIs(GI_obj) <- input |> 
    reshape2::acast(formula = as.formula(paste0(structure(GI_obj), collapse = " ~ ")), 
                    value.var = "GI")
  
  replicates(GI_obj) <- dimnames(guideGIs(GI_obj))[[which(structure(GI_obj) == "replicate")]]

  blocks(GI_obj) <- c(bio_rep = "(b\\d+)", 
            tech_rep = "(t\\d+)", 
            guide_pair = "(g\\d+)")[[block_layer]]
  
  blocks(GI_obj) <- as.character(factor(str_match(replicates(GI_obj), blocks(GI_obj))[,2]))
  ######
  GI_obj@block_description <- list(used = block_layer)
  ######
  
  return(GI_obj)
}
