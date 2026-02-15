#' @export

#create_GuideGI_arr <- function(.x, 
#                               structure, 
#                               collapse, 
#                               value_var = "GI") {
#  output <- new("GuideGI")
#  return(output)
#}

fill_guideGIs <- function(GuideGI, input, value_var = "GI") {

  #  stopifnot("Value variable not in input data" = value_var %in% colnames(.x), 
  #            "Not all dimensions given as structure found in input data!" = all(structure %in% colnames(.x)))  
  
  GuideGI@data <- reshape2::acast(
    data = input, 
    formula = as.formula(paste0(c(GuideGI@space, GuideGI@replicates), collapse = " ~ ")), 
                                 value.var = value_var)
  
  return(GuideGI)
}




GIScores <- function(input, 
                     query_col = "query_gene", 
                     lib_col = "library_gene", 
                     bio_rep_col = "bio_rep", 
                     tech_rep_col = "tech_rep", 
                     guide_col = "guide_pair", 
                     gi_col = "GI", 
                     collapse_layers = NULL, 
                     block_layer = NULL, 
                     force_fixed_pair = F, 
                     pos_agnostic = F) {
  
  #message(paste0("Collapsing layers: ", collapse_layers, collapse = ", "))

  GI_obj <- new("ScreenBase", 
                metadata = list(
                  input = input, 
                  query_col = query_col, 
                  lib_col = lib_col, 
                  bio_rep_col = bio_rep_col, 
                  tech_rep_col = tech_rep_col, 
                  guide_col = guide_col, 
                  gi_col = gi_col, 
                  force_fixed_pair = force_fixed_pair
                ))
  
  GI_obj@guideGIs@collapse <- as.character(collapse_layers)
  GI_obj@guideGIs@block_layer <- as.character(block_layer)
  
  GI_obj <- import_scores(GI_obj)
  GI_obj <- get_screen_attributes(GI_obj)
  GI_obj <- run_checks(GI_obj)
  GI_obj <- set_screenType(GI_obj)

  GI_obj@guideGIs <- GI_obj@guideGIs |> fill_guideGIs(GI_obj@metadata$input)
  GI_obj@guideGIs <- GI_obj@guideGIs |> collapse_replicates()
  GI_obj@guideGIs <- GI_obj@guideGIs |> flatten_guideGIs()
  
  print(screenReport(GI_obj))
  

  #####
  #if (is(GI_obj, "PosAgnMultiplexScreen") & pos_agnostic) {
  #  
  #  for (.r in replicates(GI_obj)) {
  #    GI_obj@guideGIs[,,.r] <- makeSymmetric(GI_obj@guideGIs[,,.r])
  #  }
  #}
  #####

  
  ##### rework?
  #GI_obj <- compute_dupCorrelation(GI_obj)
  
  #GI_obj@metadata$used_block <- block_layer
  
  return(GI_obj)
}
