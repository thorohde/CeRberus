#' @export


setGeneric("blocks", function(x) standardGeneric("blocks"))
setGeneric("blocks<-", function(x, value) standardGeneric("blocks<-"))

setGeneric("block_description", function(x) standardGeneric("block_description"))
setGeneric("block_description<-", function(x, value) standardGeneric("block_description<-"))

setGeneric("checks", function(x) standardGeneric("checks"))
setGeneric("checks<-", function(x, value) standardGeneric("checks<-"))

setGeneric("dupCorrelation", function(x) standardGeneric("dupCorrelation"))
setGeneric("dupCorrelation<-", function(x, value) standardGeneric("dupCorrelation<-"))

setGeneric("errors", function(x) standardGeneric("errors"))
setGeneric("errors<-", function(x, value) standardGeneric("errors<-"))

setGeneric("geneGIs", function(x) standardGeneric("geneGIs"))
setGeneric("geneGIs<-", function(x, value) standardGeneric("geneGIs<-"))

setGeneric("guideGIs", function(x) standardGeneric("guideGIs"))
setGeneric("guideGIs<-", function(x, value) standardGeneric("guideGIs<-"))

setGeneric("layers", function(x) standardGeneric("layers"))
setGeneric("layers<-", function(x, value) standardGeneric("layers<-"))

setGeneric("limma_models", function(x) standardGeneric("limma_models"))
setGeneric("limma_models<-", function(x, value) standardGeneric("limma_models<-"))

setGeneric("replicates", function(x) standardGeneric("replicates"))
setGeneric("replicates<-", function(x, value) standardGeneric("replicates<-"))

setGeneric("screen_attributes", function(x) standardGeneric("screen_attributes"))
setGeneric("screen_attributes<-", function(x, value) standardGeneric("screen_attributes<-"))

setGeneric("structure", function(x) standardGeneric("structure"))
setGeneric("structure<-", function(x, value) standardGeneric("structure<-"))

setGeneric("symmGeneGIs", function(x) standardGeneric("symmGeneGIs"))
setGeneric("symmGeneGIs<-", function(x, value) standardGeneric("symmGeneGIs<-"))


setGeneric("compute_dupCorrelation", function(GI_obj, ...) standardGeneric("compute_dupCorrelation"))

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