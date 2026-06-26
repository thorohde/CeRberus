# gRNA_LFC

#### Slot accessors

#slotNames("gRNA_LFC")

#### Generics for active methods

# gRNA_GIs

#### Slot accessors
# slotNames(gRNA_GI")

#### Generics for active methods

#setGeneric("fill_guideGIs", function(.x, ...) standardGeneric("fill_guideGIs"))

#' Collapse replicate dimensions in guide-level GI data
#'
#' @description
#' `collapse_replicates()` is an S4 generic used to average selected replicate
#' dimensions in guide-level genetic-interaction score containers before model
#' fitting.
#'
#' For [`gRNA_GI-class`] objects, the method reads replicate dimensions from the
#' object's `collapse` slot, verifies that all requested collapse layers are
#' valid replicate dimensions, averages over those layers with missing values
#' ignored, and updates the remaining replicate metadata.
#'
#' @param .x An object containing guide-level genetic-interaction scores.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of the same class as `.x`, with requested replicate
#'   dimensions collapsed where applicable.
#'
#' @export
setGeneric("collapse_replicates", function(.x, ...) {
  standardGeneric("collapse_replicates")
})
setGeneric("flatten_guideGIs", function(.x, ...) {
  standardGeneric("flatten_guideGIs")
})
setGeneric("compute_dupCorrelation", function(.x, ...) {
  standardGeneric("compute_dupCorrelation")
})


# ScreenBase

##### Slot accessors
# slotNames("ScreenBase")

setGeneric("checks", function(x) standardGeneric("checks"))
setGeneric("checks<-", function(x, value) standardGeneric("checks<-"))

setGeneric("dupCorrelation", function(x) standardGeneric("dupCorrelation"))
setGeneric("dupCorrelation<-", function(x, value) {
  standardGeneric("dupCorrelation<-")
})

setGeneric("errors", function(x) standardGeneric("errors"))
setGeneric("errors<-", function(x, value) standardGeneric("errors<-"))

setGeneric("geneGIs", function(x) standardGeneric("geneGIs"))
setGeneric("geneGIs<-", function(x, value) standardGeneric("geneGIs<-"))

setGeneric("guideGIs", function(x) standardGeneric("guideGIs"))
setGeneric("guideGIs<-", function(x, value) standardGeneric("guideGIs<-"))

setGeneric("guideLFCs", function(x) standardGeneric("guideLFCs"))
setGeneric("guideLFCs<-", function(x, value) standardGeneric("guideLFCs<-"))

setGeneric("limma_models", function(x) standardGeneric("limma_models"))
setGeneric("limma_models<-", function(x, value) {
  standardGeneric("limma_models<-")
})

setGeneric("screen_attr", function(x) standardGeneric("screen_attr"))
setGeneric("screen_attr<-", function(x, value) standardGeneric("screen_attr<-"))

setGeneric("symmGeneGIs", function(x) standardGeneric("symmGeneGIs"))
setGeneric("symmGeneGIs<-", function(x, value) standardGeneric("symmGeneGIs<-"))

#### Generics for active methods

#setGeneric("compute_dupCorrelation", function(GI_obj, ...) standardGeneric("compute_dupCorrelation"))
setGeneric("compute_models", function(GI_obj, ...) {
  standardGeneric("compute_models")
})
setGeneric("compute_GIs", function(GI_obj, ...) standardGeneric("compute_GIs"))
setGeneric("collect_GIs", function(GI_obj, ...) standardGeneric("collect_GIs"))
setGeneric("create_log", function(GI_obj, ...) standardGeneric("create_log"))
setGeneric("dupCorrelation_df", function(GI_obj, ...) {
  standardGeneric("dupCorrelation_df")
})
setGeneric("dupcor_df", function(GI_obj, ...) standardGeneric("dupcor_df"))
setGeneric("get_screen_attributes", function(GI_obj, ...) {
  standardGeneric("get_screen_attributes")
})
setGeneric("GI_df", function(GI_obj, ...) standardGeneric("GI_df"))
setGeneric("import_scores", function(GI_obj, ...) {
  standardGeneric("import_scores")
})
setGeneric("run_checks", function(GI_obj, ...) standardGeneric("run_checks"))
setGeneric("set_screenType", function(GI_obj, ...) {
  standardGeneric("set_screenType")
})
setGeneric("screenReport", function(GI_obj, ...) {
  standardGeneric("screenReport")
})
setGeneric("symmetry_test", function(GI_obj, ...) {
  standardGeneric("symmetry_test")
})
