#' @export


setGeneric("blocks", function(x) standardGeneric("blocks"))
setGeneric("blocks<-", function(x, value) standardGeneric("blocks<-"))

setGeneric("checks", function(x) standardGeneric("checks"))
setGeneric("checks<-", function(x, value) standardGeneric("checks<-"))

setGeneric("dupCorrelation", function(x) standardGeneric("dupCorrelation"))
setGeneric("dupCorrelation<-", function(x, value) standardGeneric("dupCorrelation<-"))

setGeneric("geneGIs", function(x) standardGeneric("geneGIs"))
setGeneric("geneGIs<-", function(x, value) standardGeneric("geneGIs<-"))

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
  "blocks", "checks", "dupCorrelation", "geneGIs", "guideGIs", 
  "replicates", "screen_attributes", "screenType", "structure"), ~ {
    
    .x2 <- paste0(.x, "<-")
    setMethod(.x, "ScreenBase", function(x) return(slot(x, .x)))
    setMethod(.x2, "ScreenBase", function(x, value) {slot(x, .x) <- value; return(x)})
    
    #  setMethod(.x, "MultiplexScreen", function(x) return(slot(x, .x)))
    #  setMethod(.x2, "MultiplexScreen", function(x, value) {slot(x, .x) <- value; return(x)})
    
    
  }
)

setGeneric("block_decision_heuristics", function(GI_obj, ...) standardGeneric("block_decision_heuristics"))
setGeneric("collapse_layer", function(GI_obj, ...) standardGeneric("collapse_layer"))
setGeneric("compute_dupcor_values", function(GI_obj, ...) standardGeneric("compute_dupcor_values"))
setGeneric("compute_GIs", function(GI_obj, ...) standardGeneric("compute_GIs"))
setGeneric("dupCorrelation_df", function(GI_obj, ...) standardGeneric("dupCorrelation_df"))
setGeneric("export_GIs", function(GI_obj, ...) standardGeneric("export_GIs"))
setGeneric("GI_df", function(GI_obj, ...) standardGeneric("GI_df"))
setGeneric("dupcor_df", function(GI_obj, ...) standardGeneric("dupcor_df"))
setGeneric("create_log", function(GI_obj, ...) standardGeneric("create_log"))

setMethod(
  "compute_dupcor_values", 
  signature = "ScreenBase", 
  function(GI_obj, sample_query = NULL) {
    
    .blocks <- blocks(GI_obj)$chosen
    if (is.null(.blocks)) {.blocks <- GI_obj@blocks$options}
    
    .blocks <- purrr::set_names(.blocks) |> purrr::map(~ GI_obj@blocks$map[,.x])
    
    if (grepl("multiplex", GI_obj@screenType)) {
      output <- if (is.null(sample_query)) {
        GI_obj@screen_attributes$query_genes
      } else {
        sample(GI_obj@screen_attributes$query_genes, min(c(sample_query, GI_obj@screen_attributes$n_query_genes)))
      }
      GI_obj@dupCorrelation <- set_names(output) |> purrr::map(~ GI_obj@guideGIs[.x,,])
    }
    
    if (grepl("fixed", GI_obj@screenType)) {
      GI_obj@dupCorrelation <- list(GI_obj@guideGIs)
    }
    
    
    
    GI_obj@dupCorrelation <- purrr::map(.blocks, \(.b) {
      GI_obj@dupCorrelation |>
        purrr::map(limma::duplicateCorrelation, block = .b) |>
        purrr::map_dbl("consensus.correlation")})
    
    return(GI_obj)
  })

setMethod(
  "collapse_layer", 
  signature = "ScreenBase", 
  function(GI_obj, .collapse) {
    
    stopifnot(.collapse %in% GI_obj@blocks$options)
    
    .old_reps <- GI_obj@blocks$options
    .new_reps <- setdiff(.old_reps, .collapse)
    
    output <- input |> flatten_array(structure(GI_obj), "GI")
    
    output[, c(.old_reps) := data.table::tstrsplit(replicate, "_")]
    
    output[, replicate := do.call(paste, c(data.table::.SD, sep = "_")), .SDcols = .new_reps]
    
    output <- output |> 
      reshape2::acast(formula = as.formula(paste0(structure(GI_obj), collapse = " ~ ")), 
                      value.var = "GI", 
                      fun.aggregate = mean)
    
    GI_obj@blocks$collapsed <- c(GI_obj@blocks$collapsed, .collapse)
    
    replicates(GI_obj) <- dimnames(guideGIs(GI_obj))[[which(structure(GI_obj) == "replicate")]]
    
    blocks(GI_obj) <- list(map = map_replicate_layers(replicates(GI_obj)), 
                           all = colnames(map_replicate_layers(replicates(GI_obj))), 
                           options = colnames(map_replicate_layers(replicates(GI_obj))), 
                           collapsed = NULL, 
                           chosen = NULL)
    return(GI_obj)
  })

setMethod(
  "block_decision_heuristics", 
  signature = "ScreenBase", 
  function(GI_obj) {
    
    while(is.null(GI_obj@blocks$chosen) & 
          ncol(GI_obj@blocks$map) > 1) {
      
      # This should find the block with the highest in-block correlation. 
      # If no positive block is found, the lowest correlating layer gets collapsed, 
      # and the process is repeated. 
      
      .dc <- GI_obj@dupCorrelation
      
      if (!all(colnames(GI_obj@blocks$map) %in% names(.dc))) {
        
        GI_obj <- compute_dupcor_values(GI_obj, sample_query = 20)
        .dc <- GI_obj@dupCorrelation
      }
      
      .dc <- purrr::map_dbl(.dc, median, na.rm = T)
      
      .highest <- .dc[which.max(.dc)]
      
      if (.highest >= 0) {
        GI_obj@blocks$chosen <- names(.highest)
      } else {
        GI_obj <- collapse_layer(GI_obj, names(which.min(.dc)))
      }}
    return(GI_obj)
  }
)

setMethod(
  "compute_GIs", 
  signature = "ScreenBase", 
  function(GI_obj, FDR_method) {
    
    stopifnot("FDR method required." = !missing(FDR_method))
    .block <- GI_obj@blocks$chosen
    
    if (grepl("multiplex", GI_obj@screenType)) {
      
      geneGIs(GI_obj) <- set_names(GI_obj@screen_attributes$query_genes) |> 
        map(purrr::safely(\(.g) {
          .fit <- limma::lmFit(object = GI_obj@guideGIs[.g,,], 
                               block = GI_obj@blocks$map[, GI_obj@blocks$chosen], 
                               correlation = GI_obj@dupCorrelation[[GI_obj@blocks$chosen]][[.g]])
          
          .efit <- limma::eBayes(.fit)
          
          data.frame(query_gene = rep(.g, GI_obj@screen_attributes$n_lib_genes), 
                     library_gene = GI_obj@screen_attributes$library_genes, 
                     GI = .fit$Amean, 
                     pval = .efit$p.value[, 1], 
                     FDR = stats::p.adjust(.efit$p.value[, 1], method = FDR_method))
        })) |> purrr::map("result") |> 
        data.table::rbindlist(fill = T) |>
        data.table::melt.data.table(measure.vars = c("GI", "pval", "FDR")) |>
        reshape2::acast(formula = as.formula("query_gene ~ library_gene ~ variable"), 
                        value.var = "value")
      
    } else {
      geneGIs(GI_obj) <- {
        .fit <- limma::lmFit(object = GI_obj@guideGIs, 
                             block = GI_obj@blocks$map[, GI_obj@blocks$chosen], 
                             correlation = GI_obj@dupCorrelation[[GI_obj@blocks$chosen]])
        
        
        .efit <- limma::eBayes(.fit)
        
        output <- list(rownames(GI_obj@guideGIs), 
                       c("GI", "pval", "FDR"))
        
        output <- array(data = NA, 
                        dim = map_int(output, length), 
                        dimnames = output)
        
        output[,"GI"] <- .fit$Amean
        output[,"pval"] <- .efit$p.value[, 1]
        output[,"FDR"] <- stats::p.adjust(.efit$p.value[, 1], method = FDR_method)
        
        output
      }}
    
    return(GI_obj)
  })

setMethod(
  "GI_df", 
  signature = "ScreenBase", 
  function(GI_obj) {
    
    if (grepl("multiplex", GI_obj@screenType)) {
      output <- GI_obj@geneGIs |> 
        flatten_array(dnames = c("query_gene", "library_gene", "variable"), 
                      value.name = "value") |>
        data.table::dcast(formula = query_gene + library_gene ~ variable, 
                          value.var = "value")
      output <- output[, .SD, .SDcols = c("query_gene", "library_gene", "GI", "pval", "FDR")]
    }
    
    if (grepl("fixed", GI_obj@screenType)) {
      output <- data.table::data.table(pair = rownames(GI_obj@geneGIs), GI_obj@geneGIs)
    }
    
    return(output)
    
  })


setMethod(
  "dupCorrelation_df", 
  signature = "ScreenBase", 
  function(GI_obj) {
    if (grepl("multiplex", GI_obj@screenType)) {
      output <- data.table::data.table(
        query_gene = GI_obj@screen_attributes$query_genes, 
        data.frame(GI_obj@dupCorrelation))
    }
    
    if (grepl("fixed", GI_obj@screenType)) {
      output <- data.table(data.frame(GI_obj@dupCorrelation))
    }
    
    
    return(output)
    
  })

setMethod(
  "create_log", 
  signature = "ScreenBase", 
  function(GI_obj) {
    list(stringr::str_glue("CRISPR screen with {GI_obj@screen_attributes$n_query_genes} query genes and {GI_obj@screen_attributes$n_lib_genes} library genes."), 
         stringr::str_glue("Identified screen design: {GI_obj@screenType}."), 
         "", 
         stringr::str_glue("Possible replicate layers: {paste(GI_obj@blocks$options, collapse = ', ')}"), 
         stringr::str_glue("Possible layer used for p-values: {GI_obj@blocks$chosen}"), 
         stringr::str_glue("Median duplicate correlation levels: "), 
         "", 
         stringr::str_glue("... and many other useful log information.")
    ) |> paste(collapse = "\n")
    
  })

