#' @export

setClass("GuideGI", 
         slots = list(
           "data" = "array", 
           "space" = "character", 
           "replicates" = "character", 
           "block_layer" = "character", 
           "blocks" = "character", 
           "use_blocks" = "logical", 
           "block_description" = "character", 
           "collapse" = "character")
         )

setClass("ScreenBase", 
         slots = list(
           "guideGIs" = "GuideGI", 
           "limma_models" = "list", 
           "geneGIs" = "array", 
           "screen_attr" = "list", 
           #"blocks" = "character", 
           #"block_description" = "list", 
           "dupCorrelation" = "numeric",
           "metadata" = "list", 
           "checks" = "list",
           "errors" = "list"))

setClass("FixedPairScreen", contains = "ScreenBase")
setClass("MultiplexScreen", contains = "ScreenBase")
setClass("PosAgnMultiplexScreen", contains = "MultiplexScreen", slots = list("symmGeneGIs" = "data.table"))


# Need to be called without a loop to work

setAs(from = "ScreenBase", to = "FixedPairScreen", 
      function(from) {
        obj <- new("FixedPairScreen")
        for (s in slotNames("ScreenBase")) {slot(obj, s) <- slot(from, s)}
        class(obj) <- "FixedPairScreen"
        return(obj)})

setAs(from = "ScreenBase", to = "MultiplexScreen", 
      function(from) {
        obj <- new("MultiplexScreen")
        for (s in slotNames("ScreenBase")) {slot(obj, s) <- slot(from, s)}
        class(obj) <- "MultiplexScreen"
        return(obj)})

setAs(from = "ScreenBase", to = "PosAgnMultiplexScreen", 
      function(from) {
        obj <- new("PosAgnMultiplexScreen")
        for (s in slotNames("ScreenBase")) {slot(obj, s) <- slot(from, s)}
        class(obj) <- "PosAgnMultiplexScreen"
        return(obj)})