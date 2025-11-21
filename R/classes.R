#' @export

setClass("ScreenBase", 
         slots = c(
           "screenType" = "character", 
           "guideGIs" = "array", 
           #"guideGIsSymmetric" = "array", 
           "geneGIs" = "array", 
           "geneGIsSymmetric" = "array", 
           "structure" = "character", 
           "replicates" = "character",  
           "screen_attributes" = "list", 
           "blocks" = "character", 
           "block_description" = "list", 
           "dupCorrelation" = "numeric",
           "checks" = "list",
           "errors" = "list")
)

#setClass("FixedPairScreen", contains = "ScreenBase")
