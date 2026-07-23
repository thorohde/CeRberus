# gRNA_LFC

#### Slot accessors

#### Generics for active methods

# gRNA_GIs

#### Slot accessors

#### Generics for active methods

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
setGeneric("flatten_guide_gis", function(.x, ...) {
  standardGeneric("flatten_guide_gis")
})
setGeneric("compute_dup_correlation", function(.x, ...) {
  standardGeneric("compute_dup_correlation")
})


# ScreenBase

##### Slot accessors

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

setGeneric("compute_models", function(gi_obj, ...) {
  standardGeneric("compute_models")
})
setGeneric("collect_gis", function(gi_obj, ...) standardGeneric("collect_gis"))
setGeneric("create_log", function(gi_obj, ...) standardGeneric("create_log"))
setGeneric("dup_correlation_df", function(gi_obj, ...) {
  standardGeneric("dup_correlation_df")
})
setGeneric("get_screen_attributes", function(gi_obj, ...) {
  standardGeneric("get_screen_attributes")
})
setGeneric("gi_df", function(gi_obj, ...) standardGeneric("gi_df"))
setGeneric("import_scores", function(gi_obj, ...) {
  standardGeneric("import_scores")
})
setGeneric("run_checks", function(gi_obj, ...) standardGeneric("run_checks"))
setGeneric("set_screen_type", function(gi_obj, ...) {
  standardGeneric("set_screen_type")
})
setGeneric("screen_report", function(gi_obj, ...) {
  standardGeneric("screen_report")
})
setGeneric("symmetry_test", function(gi_obj, ...) {
  standardGeneric("symmetry_test")
})
