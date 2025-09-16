#' @export

setClass("GIScores", 
         slots = c("guideGIs" = "array", 
                   "dupcor" = "list",
                   "layers" = "matrix", 
                   "dim_description" = "integer", 
                   "screen_attributes" = "list", 
                   "checks" = "list", 
                   "screen_type" = "character")
)


