#' @include generics.R
#' @export

setMethod("dim_description", "GIScores", function(x) return(x@dim_description))
setMethod("dim_description<-", "GIScores", function(x, value) {x@dim_description <- value; return(x)})

setMethod("guideGIs", "GIScores", function(x) return(x@guideGIs))
setMethod("guideGIs<-", "GIScores", function(x, value) {x@guideGIs <- value; return(x)})

setMethod("layers", "GIScores", function(x) return(x@layers))
setMethod("layers<-", "GIScores", function(x, value) {x@layers <- value; return(x)})

setMethod("screen_attributes", "GIScores", function(x) return(x@screen_attributes))
setMethod("screen_attributes<-", "GIScores", function(x, value) {x@screen_attributes <- value; return(x)})


setMethod("checks", "GIScores", function(x) return(x@checks))
setMethod("checks<-", "GIScores", function(x, value) {x@checks <- value; return(x)})


setMethod("screen_type", "GIScores", function(x) return(x@screen_type))
setMethod("screen_type<-", "GIScores", function(x, value) {x@screen_type <- value; return(x)})
