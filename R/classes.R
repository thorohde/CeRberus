#' @export

setClass("ScreenBase", 
         slots = c(
           "screenType" = "character", 
           "guideGIs" = "array", 
           "geneGIs" = "array", 
           "geneGIsSymmetric" = "list", 
           "structure" = "character", 
           "replicates" = "character",  
           "screen_attributes" = "list", 
           "blocks" = "character", 
           "block_description" = "list", 
           "dupCorrelation" = "numeric",
           "checks" = "list"
           )
)

#setClass("FixedPairScreen", contains = "ScreenBase")
