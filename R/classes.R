#' @export

setClass("ScreenBase", 
         slots = c(
           "screenType" = "character", 
           "guideGIs" = "array", 
           "geneGIs" = "list", 
           "geneGIsSymmetric" = "list", 
           "structure" = "character", 
           "replicates" = "character",  
           "screen_attributes" = "list", 
           "blocks" = "list", 
           "dupCorrelation" = "list",
           "checks" = "list"
           )
)

#setClass("FixedPairScreen", contains = "ScreenBase")
