#' @export


# Slot accessors for guideGIs:


##### Active methods for guideGIs:
#setGeneric("fill_guideGIs", function(.x, ...) standardGeneric("fill_guideGIs"))
setGeneric("collapse_replicates", function(.x, ...) standardGeneric("collapse_replicates"))
setGeneric("flatten_guideGIs", function(.x, ...) standardGeneric("flatten_guideGIs"))
setGeneric("compute_dupCorrelation", function(.x, ...) standardGeneric("compute_dupCorrelation"))




# Slot accessors for GI_objects:

walk(c(#"blocks", 
  #"block_description", 
  "checks", 
  "dupCorrelation", 
  "errors", 
  "geneGIs", 
  "guideGIs", 
  "limma_models", 
  #"replicates", 
  "screen_attr", 
  "symmGeneGIs"
), function(name) {
  setGeneric(name, eval(substitute(function(x) standardGeneric(NAME), list(NAME = name))))
  setGeneric(paste0(name, "<-"), eval(substitute(function(x, value) standardGeneric(NAME),
                                                 list(NAME = paste0(name, "<-")))))})





##### Active methods for GI_objects:
#setGeneric("compute_dupCorrelation", function(GI_obj, ...) standardGeneric("compute_dupCorrelation"))
setGeneric("compute_models", function(GI_obj, ...) standardGeneric("compute_models"))
setGeneric("compute_GIs", function(GI_obj, ...) standardGeneric("compute_GIs"))
setGeneric("collect_GIs", function(GI_obj, ...) standardGeneric("collect_GIs"))
setGeneric("create_log", function(GI_obj, ...) standardGeneric("create_log"))
setGeneric("dupCorrelation_df", function(GI_obj, ...) standardGeneric("dupCorrelation_df"))
setGeneric("dupcor_df", function(GI_obj, ...) standardGeneric("dupcor_df"))
setGeneric("get_screen_attributes", function(GI_obj, ...) standardGeneric("get_screen_attributes"))
setGeneric("GI_df", function(GI_obj, ...) standardGeneric("GI_df"))
setGeneric("import_scores", function(GI_obj, ...) standardGeneric("import_scores"))
setGeneric("run_checks", function(GI_obj, ...) standardGeneric("run_checks"))
setGeneric("set_screenType", function(GI_obj, ...) standardGeneric("set_screenType"))
setGeneric("screenReport", function(GI_obj, ...) standardGeneric("screenReport"))
setGeneric("symmetry_test", function(GI_obj, ...) standardGeneric("symmetry_test"))
