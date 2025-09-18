#' @export

setClass("ScreenBase", 
         slots = c(
           "screenType" = "character", 
           "guideGIs" = "array", 
           "geneGIs" = "array", 
           "structure" = "character", 
           "replicates" = "character",  
           "screen_attributes" = "list", 
           "blocks" = "list", 
           "dupCorrelation" = "list",
           "checks" = "list"
           )
)

#setClass("FixedPairScreen", contains = "ScreenBase")
#setClass("MultiplexScreen", contains = "ScreenBase")
