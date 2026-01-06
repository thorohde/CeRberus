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

setGeneric("geneGIsSymmetric", function(x) standardGeneric("geneGIsSymmetric"))
setGeneric("geneGIsSymmetric<-", function(x, value) standardGeneric("geneGIsSymmetric<-"))

setGeneric("guideGIs", function(x) standardGeneric("guideGIs"))
setGeneric("guideGIs<-", function(x, value) standardGeneric("guideGIs<-"))

setGeneric("layers", function(x) standardGeneric("layers"))
setGeneric("layers<-", function(x, value) standardGeneric("layers<-"))

setGeneric("replicates", function(x) standardGeneric("replicates"))
setGeneric("replicates<-", function(x, value) standardGeneric("replicates<-"))

setGeneric("screen_attributes", function(x) standardGeneric("screen_attributes"))
setGeneric("screen_attributes<-", function(x, value) standardGeneric("screen_attributes<-"))

setGeneric("screenType", function(x) standardGeneric("screenType"))
setGeneric("screenType<-", function(x, value) standardGeneric("screenType<-"))

setGeneric("structure", function(x) standardGeneric("structure"))
setGeneric("structure<-", function(x, value) standardGeneric("structure<-"))





purrr::walk(c(
  "blocks", "block_description", "checks", "dupCorrelation", 
  "errors", "geneGIs", "geneGIsSymmetric", "guideGIs", 
  "replicates", "screen_attributes", "screenType", "structure"), ~ {
    
    .x2 <- paste0(.x, "<-")
    setMethod(.x, "ScreenBase", function(x) return(slot(x, .x)))
    setMethod(.x2, "ScreenBase", function(x, value) {slot(x, .x) <- value; return(x)})
    #  setMethod(.x, "MultiplexScreen", function(x) return(slot(x, .x)))
    #  setMethod(.x2, "MultiplexScreen", function(x, value) {slot(x, .x) <- value; return(x)})
    
  }
)

setGeneric("compute_dupCorrelation", function(GI_obj, ...) standardGeneric("compute_dupCorrelation"))
setGeneric("compute_GIs", function(GI_obj, ...) standardGeneric("compute_GIs"))
setGeneric("compute_symmetric_GIs", function(GI_obj, ...) standardGeneric("compute_symmetric_GIs"))
setGeneric("dupCorrelation_df", function(GI_obj, ...) standardGeneric("dupCorrelation_df"))
setGeneric("GI_df", function(GI_obj, ...) standardGeneric("GI_df"))
setGeneric("symmetricGI_df", function(GI_obj, ...) standardGeneric("symmetricGI_df"))
setGeneric("dupcor_df", function(GI_obj, ...) standardGeneric("dupcor_df"))
setGeneric("create_log", function(GI_obj, ...) standardGeneric("create_log"))



setMethod("compute_dupCorrelation", signature = "ScreenBase", 
          function(GI_obj) {
            if (grepl("fixed", GI_obj@screenType)) {
              
              dupCorrelation(GI_obj) <- suppressWarnings(limma::duplicateCorrelation(
                object = GI_obj@guideGIs, 
                block = blocks(GI_obj)))$consensus.correlation
            }
            
            if (grepl("multiplex", GI_obj@screenType)) {
              
              output <- purrr::set_names(GI_obj@screen_attributes$query_genes)
              
              output <- output |> purrr::map(\(.g) GI_obj@guideGIs[.g,,])
              
              output <- output |>
                purrr::map(\(.g) suppressWarnings(limma::duplicateCorrelation(object = .g, 
                                                                              block = blocks(GI_obj)))) |>
                purrr::map_dbl("consensus.correlation")
              
              GI_obj@dupCorrelation <- output
            }
            return(GI_obj)
          })




setMethod(
  "compute_GIs", 
  signature = "ScreenBase", 
  function(GI_obj, FDR_method = "BH") {
    
    stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)
    
    #stopifnot("DupCorrelation values required for all blocks!" = all(.blocks %in% names(GI_obj@dupCorrelation)))
    
    if (grepl("multiplex", GI_obj@screenType)) {
      
      output <- set_names(GI_obj@screen_attributes$query_genes)
      
      output <- output |> 
        map(safely(\(.g) {
          
          .obj <- GI_obj@guideGIs[.g,,]
          
          .dcor <- dupCorrelation(GI_obj)[[.g]]
          
          .fit <- limma::lmFit(object = .obj, 
                               block = blocks(GI_obj), 
                               correlation = .dcor)
          
          .efit <- limma::eBayes(.fit)
          
          .obj <- data.frame(library_gene = GI_obj@screen_attributes$library_genes, 
                             GI = .fit$Amean, 
                             pval = .efit$p.value[, 1], 
                             FDR = stats::p.adjust(.efit$p.value[, 1], method = FDR_method))
          
          return(.obj)
        }))
      
      errors(GI_obj)$query_genes_not_usable <- map(output, "result") |> keep(is.null) |> names()
      errors(GI_obj)$GI_computation_errors <- map(output, "error")
      
      #  map("error") |> 
      #  compact() |> 
      #  map_chr(~ .x$message)
      
      if (length(errors(GI_obj)$query_genes_not_usable) > 0) {
        warning(str_c("Failed computing GIs for ", length(errors(GI_obj)$query_genes_not_usable), " genes."))
      }
      
      output <- output |> 
        map("result") |>
        modify_if(is.null, 
                  \(.obj) {data.frame(library_gene = GI_obj@screen_attributes$library_genes, 
                                      GI = NA, pval = NA, FDR = NA)})
      
      
      
      
      
      
      
      output <- output |> 
        imap(~ {.x$query_gene <- .y; .x}) |>
        data.table::rbindlist(fill = T) |>
        data.table::melt.data.table(measure.vars = c("GI", "pval", "FDR")) |>
        reshape2::acast(formula = as.formula("query_gene ~ library_gene ~ variable"), 
                        value.var = "value", drop = F)
      
    } else {
      .fit <- limma::lmFit(object = GI_obj@guideGIs, 
                           block = blocks(GI_obj), 
                           correlation = dupCorrelation(GI_obj))
      
      .efit <- limma::eBayes(.fit)
      
      output <- list(rownames(GI_obj@guideGIs), 
                     c("GI", "pval", "FDR"))
      
      output <- array(data = NA, 
                      dim = purrr::map_int(output, length), 
                      dimnames = output)
      
      output[,"GI"] <- .fit$Amean
      output[,"pval"] <- .efit$p.value[, 1]
      output[,"FDR"] <- stats::p.adjust(.efit$p.value[, 1], method = FDR_method)
      
    }
    
    geneGIs(GI_obj) <- output
    
    if (length(setdiff(rownames(guideGIs(GI_obj)), rownames(geneGIs(GI_obj)))) != 0 | 
        length(setdiff(colnames(guideGIs(GI_obj)), colnames(geneGIs(GI_obj)))) != 0) {
      
      warning("Some genes were lost!")
      
      print(str(guideGIs(GI_obj)))
      print(str(geneGIs(GI_obj)))
      
      print(setdiff(rownames(guideGIs(GI_obj)), rownames(geneGIs(GI_obj))))
      print(setdiff(colnames(guideGIs(GI_obj)), colnames(geneGIs(GI_obj))))
      
    }
    
    
    
    return(GI_obj)
  })

setMethod("compute_symmetric_GIs", 
          signature = "ScreenBase", 
          function(GI_obj, FDR_method = "BH") {
            
            if (!grepl("position.agnostic", screenType(GI_obj))) {
              warning("The given screen is not position-agnostic.")
              return(GI_obj)
            }
            
            .x <- list(screen_attributes(GI_obj)$unique_pairs, c("GI", "GI_z", "pval", "FDR"))
            .x <- array(data = NA, 
                        dim = purrr::map_int(.x, length), 
                        dimnames = .x)
            
            .arr <- geneGIs(GI_obj)
            
            .x[,"GI"] <- gather_symmetric_scores(pairs = rownames(.x), 
                                                 .arr = .arr[,,"GI"])
            .x[,"GI_z"] <- z_transform(.x[,"GI"])
            .x[,"pval"] <- gather_symmetric_scores(pairs = rownames(.x), 
                                                   .arr = .arr[,,"pval"])
            
            .x[,"FDR"] <- balanced_FDR(pairs = rownames(.x), 
                                       pval_array = .arr[,,"pval"], 
                                       fdr_method = FDR_method)
            
            geneGIsSymmetric(GI_obj) <- .x
            
            return(GI_obj)
          })

setMethod(
  "GI_df", 
  signature = "ScreenBase", 
  function(GI_obj) {
    
    if (grepl("multiplex", GI_obj@screenType)) {
      output <- GI_obj@geneGIs
      
      output <- output |> 
        flatten_array(dnames = c("query_gene", "library_gene", "variable"), 
                      value.name = "value") |>
        data.table::dcast(formula = query_gene + library_gene ~ variable, 
                          value.var = "value")
      output <- output[, .SD, .SDcols = c("query_gene", "library_gene", "GI", "pval", "FDR")]
    }
    
    if (grepl("fixed", GI_obj@screenType)) {
      output <- GI_obj@geneGIs
      
      output <- data.table::data.table(pair = rownames(output), 
                                       query_gene = str_split_i(rownames(output), ";", 1), 
                                       library_gene = str_split_i(rownames(output), ";", 2), 
                                       output)
    }
    return(output)
  })



setMethod(
  "symmetricGI_df", 
  signature = "ScreenBase", 
  function(GI_obj) {
    
    if (GI_obj@screenType == "multiplex.symmetric.position.agnostic") {
      
      
      output <- GI_obj@geneGIsSymmetric
      
      output <- data.table::data.table(pair = rownames(output), 
                                       query_gene = str_split_i(rownames(output), ";", 1), 
                                       library_gene = str_split_i(rownames(output), ";", 2), 
                                       output)
      
      return(output)
    } else {
      return(NULL)
    }
  })






setMethod(
  "dupCorrelation_df", 
  signature = "ScreenBase", 
  function(GI_obj) {
    output <- data.table(data.frame(dupcor = GI_obj@dupCorrelation))
    return(output)
  })

setMethod(
  "create_log", 
  signature = "ScreenBase", 
  function(GI_obj) {
    list(stringr::str_glue("CRISPR screen with {GI_obj@screen_attributes$n_query_genes} query genes and {GI_obj@screen_attributes$n_lib_genes} library genes."), 
         stringr::str_glue("Identified screen design: {GI_obj@screenType}."), 
         "", 
         #stringr::str_glue("Possible replicate layers: {paste(GI_obj@blocks$options, collapse = ', ')}"), 
         #stringr::str_glue("Possible layer used for p-values: {GI_obj@blocks$chosen}"), 
         stringr::str_glue("Median duplicate correlation levels: "), 
         "", 
         stringr::str_glue("... and many other useful log information.")
    ) |> paste(collapse = "\n")
    
  })



