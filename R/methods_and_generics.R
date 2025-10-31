#' @export


setGeneric("blocks", function(x) standardGeneric("blocks"))
setGeneric("blocks<-", function(x, value) standardGeneric("blocks<-"))

setGeneric("block_description", function(x) standardGeneric("block_description"))
setGeneric("block_description<-", function(x, value) standardGeneric("block_description<-"))

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
  "blocks", "block_description", "checks", "dupCorrelation", "geneGIs", "guideGIs", 
  "replicates", "screen_attributes", "screenType", "structure"), ~ {
    
    .x2 <- paste0(.x, "<-")
    setMethod(.x, "ScreenBase", function(x) return(slot(x, .x)))
    setMethod(.x2, "ScreenBase", function(x, value) {slot(x, .x) <- value; return(x)})
    
    #  setMethod(.x, "MultiplexScreen", function(x) return(slot(x, .x)))
    #  setMethod(.x2, "MultiplexScreen", function(x, value) {slot(x, .x) <- value; return(x)})
    
    
  }
)

#setGeneric("add_collapsed_layers", function(GI_obj, ...) standardGeneric("add_collapsed_layers"))
setGeneric("block_decision_heuristics", function(GI_obj, ...) standardGeneric("block_decision_heuristics"))
#setGeneric("collapse_layer", function(GI_obj, ...) standardGeneric("collapse_layer"))
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
    
    if (grepl("multiplex", GI_obj@screenType)) {
      output <- set_names(GI_obj@screen_attributes$query_genes)
      
      if (!is.null(sample_query)) {
        output <- sample(output, min(c(sample_query, GI_obj@screen_attributes$n_query_genes)))
      }
      output <- output |> purrr::map(\(.g) GI_obj@guideGIs[.g,,])
    }
    
    if (grepl("fixed", GI_obj@screenType)) {
      output <- list(GI_obj@guideGIs)
    }
    
    output <- output |>
      purrr::map(\(.g) suppressWarnings(limma::duplicateCorrelation(object = .g, 
                                                                block = blocks(GI_obj)))) |>
      purrr::map_dbl("consensus.correlation")
    
    GI_obj@dupCorrelation <- output
    
    return(GI_obj)
  })


## rework!
if (F) {
  setMethod(
  "block_decision_heuristics", 
  signature = "ScreenBase", 
  function(GI_obj) {
    
    while(is.null(GI_obj@blocks$chosen) & 
          length(GI_obj@blocks$map) > 1) {
      
      # This should find the block with the highest in-block correlation. 
      # If no positive block is found, the lowest correlating layer gets collapsed, 
      # and the process is repeated. 
      
      .dc <- GI_obj@dupCorrelation
      
      if (!all(names(GI_obj@blocks$map) %in% names(.dc))) {
        
        GI_obj <- compute_dupcor_values(GI_obj, sample_query = 20)
        #any(dupCorrelation(.d) |> map(mean, na.rm = T) > 0.05)
        
        .dc <- GI_obj@dupCorrelation
      }
      
      .dc <- purrr::map_dbl(.dc, median, na.rm = T)
      
      .highest <- .dc[which.max(.dc)]
      
      if (.highest >= 0) {
        GI_obj@blocks$chosen <- names(.highest)
      } else {
        print(stringr::str_glue("Collapsing {names(which.min(.dc))}"))
        #GI_obj <- collapse_layer(GI_obj, names(which.min(.dc)))
      }}
    return(GI_obj)
  }
)
}
setMethod(
  "compute_GIs", 
  signature = "ScreenBase", 
  function(GI_obj, FDR_method = "BH") {
    
    stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)
    
    #stopifnot("DupCorrelation values required for all blocks!" = all(.blocks %in% names(GI_obj@dupCorrelation)))
    
    if (grepl("multiplex", GI_obj@screenType)) {
      
      
      output <- set_names(GI_obj@screen_attributes$query_genes) |> 
        map(purrr::safely(\(.g) {
          .fit <- limma::lmFit(object = GI_obj@guideGIs[.g,,], 
                               block = blocks(GI_obj), 
                               correlation = dupCorrelation(GI_obj)[[.g]])
          
          .efit <- limma::eBayes(.fit)
          
          data.frame(query_gene = .g, 
                     library_gene = rownames(.fit), 
                     GI = .fit$Amean, 
                     pval = .efit$p.value[, 1], 
                     FDR = stats::p.adjust(.efit$p.value[, 1], method = FDR_method))
        })) |> 
        purrr::map("result") |> 
        data.table::rbindlist(fill = T) |>
        data.table::melt.data.table(measure.vars = c("GI", "pval", "FDR")) |>
        reshape2::acast(formula = as.formula("query_gene ~ library_gene ~ variable"), 
                        value.var = "value")
      
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
    
    return(GI_obj)
  })

# rework
compute_symmetric_GIs <- function(GI_arr, FDR_method = "BH") {
  
  stopifnot(
    "GI Array is not symmetric!" = 
      cor(x = as.vector(GI_arr[,,"GI"]), 
          y = as.vector(t(GI_arr[,,"GI"])), 
          use = "pairwise.complete.obs") >= 0.98, 
    
    "pval array is not symmetric!" = 
      cor(x = as.vector(GI_arr[,,"pval"]), 
          y = as.vector(t(GI_arr[,,"pval"])), 
          use = "pairwise.complete.obs") >= 0.9)
  
  gather_symmetric_scores <- function(pairs, .arr, sep = ";") {
    if (!isSymmetric(.arr)) {
      warning(stringr::str_c("Input array is asymmetric! ABBA: ", 
                             round(cor(x = as.vector(.arr), 
                                       y = as.vector(t(.arr)), 
                                       use = "pairwise.complete.obs"), 3)))}
    
    genes1 <- str_split_i(pairs, sep, 1)
    genes2 <- str_split_i(pairs, sep, 2)
    return(purrr::map2_dbl(genes1, genes2, \(.g1, .g2) {.arr[.g1,.g2]}))}
  
  balanced_FDR <- \(pairs, pval_array, fdr_method) {
    
    pair_template <- \(.pair) {
      .g1 <- str_split_i(.pair, ";", 1)
      .g2 <- str_split_i(.pair, ";", 2)
      .d <- data.table(pair = .pair, 
                       g1 = c(rep(.g1, ncol(pval_array)), setdiff(rownames(pval_array), .g1)), 
                       g2 = c(colnames(pval_array), rep(.g2, nrow(pval_array)-1)))
      .d[, pair2 := stringr::str_c(g1, ";", g2)]
      .d
    }
    
    return(purrr::map_dbl(pairs, \(pair) {
      .d <- pair_template(pair)
      .d[, pval := purrr::map2_dbl(g1, g2, ~ pval_array[.x,.y])]
      .d[, fdr := stats::p.adjust(pval, method = fdr_method)]
      .d[pair == pair2, get("fdr")]}))
  }
  
  
  all_pairs <- data.table::CJ(g1 = rownames(GI_arr), g2 = colnames(GI_arr))
  all_pairs <- all_pairs[, `:=`(pair = stringr::str_glue("{g1};{g2}"), 
                                sorted_pair = sort_gene_pairs(g1, g2))]
  
  .x <- list(all_pairs[, unique(sorted_pair)], 
             c("GI", "GI_z", "pval", "fdr"))
  
  .x <- array(data = NA, 
              dim = purrr::map_int(.x, length), 
              dimnames = .x)
  
  .x[,"GI"] <- gather_symmetric_scores(pairs = rownames(.x), 
                                       .arr = GI_arr[,,"GI"])
  .x[,"GI_z"] <- z_transform(.x[,"GI"])
  .x[,"pval"] <- gather_symmetric_scores(pairs = rownames(.x), 
                                         .arr = GI_arr[,,"pval"])
  
  .x[,"fdr"] <- balanced_FDR(pairs = rownames(.x), 
                             pval_array = GI_arr[,,"pval"], 
                             fdr_method = FDR_method)
  
  return(.x)
}





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
      output <- data.table::data.table(pair = rownames(GI_obj@geneGIs), 
                                       query_gene = str_split_i(rownames(GI_obj@geneGIs), ";", 1), 
                                       library_gene = str_split_i(rownames(GI_obj@geneGIs), ";", 2), 
                                       GI_obj@geneGIs)
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
        data.frame(dupcor = GI_obj@dupCorrelation))
    }
    if (grepl("fixed", GI_obj@screenType)) {
      output <- data.table(data.frame(dupcor = GI_obj@dupCorrelation))
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
         #stringr::str_glue("Possible replicate layers: {paste(GI_obj@blocks$options, collapse = ', ')}"), 
         #stringr::str_glue("Possible layer used for p-values: {GI_obj@blocks$chosen}"), 
         stringr::str_glue("Median duplicate correlation levels: "), 
         "", 
         stringr::str_glue("... and many other useful log information.")
    ) |> paste(collapse = "\n")
    
  })

