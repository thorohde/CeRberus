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
                     force_fixed_pair = F, 
                     pos_agnostic = F#, 
                     #individual_guide_dupcor = T
) {
  
  #message(paste0("Collapsing layers: ", collapse_layers, collapse = ", "))
  #message(paste0(block_layer, " used as blocks"))
  
  GI_obj <- new("ScreenBase", 
                metadata = list(
                  input = input, 
                  query_col = query_col, 
                  lib_col = lib_col, 
                  bio_rep_col = bio_rep_col, 
                  tech_rep_col = tech_rep_col, 
                  guide_col = guide_col, 
                  gi_col = gi_col, 
                  collapse_layers = collapse_layers, 
                  force_fixed_pair = force_fixed_pair
                ))
  
  GI_obj <- import_scores(GI_obj)
  
  GI_obj <- get_screen_attributes(GI_obj)
  
  GI_obj <- run_checks(GI_obj)
  
  GI_obj <- set_screenType(GI_obj)
  
  print(screenReport(GI_obj))
  
    
  if (is(GI_obj, "PosAgnMultiplexScreen") & pos_agnostic) {
    
    for (.r in replicates(GI_obj)) {
      guideGIs(GI_obj)[,,.r] <- makeSymmetric(guideGIs(GI_obj)[,,.r])
    }
  }
  
  if (!is.null(block_layer)) {
    blocks(GI_obj) <- c(bio_rep = "(b\\d+)", 
                        tech_rep = "(t\\d+)", 
                        guide_pair = "(g\\d+)")[[block_layer]]
    
    blocks(GI_obj) <- as.character(factor(str_match(replicates(GI_obj), blocks(GI_obj))[,2]))
  } else {
    warning("Not defining a block structure is not recommended.")
  }
  
  ##### rework?
  GI_obj <- compute_dupCorrelation(GI_obj)
  
  GI_obj@block_description <- list(used = block_layer)
  
  return(GI_obj)
}
