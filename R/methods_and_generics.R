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

setGeneric("add_collapsed_layers", function(GI_obj, ...) standardGeneric("add_collapsed_layers"))
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
  function(GI_obj, .blocks = GI_obj@blocks$use, sample_query = NULL) {
    
    if (is.null(.blocks)) {.blocks <- GI_obj@blocks$all}
    
    .blocks <- GI_obj@blocks$map[.blocks]
    
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
        purrr::map(~ suppressWarnings(limma::duplicateCorrelation(object = .x[,names(.b)], block = .b))) |>
        purrr::map_dbl("consensus.correlation")})
    
    
    return(GI_obj)
  })

#setMethod(
#    "collapse_layer", 
#    signature = "ScreenBase", 
#    function(GI_obj, .collapse) {

#      stopifnot(.collapse %in% GI_obj@blocks$options)
#      
#      .old_reps <- GI_obj@blocks$options
#      .new_reps <- setdiff(.old_reps, .collapse)
#      
#      output <- GI_obj@guideGIs |> flatten_array(structure(GI_obj), "GI")
#      
#      output[, c(.old_reps) := data.table::tstrsplit(replicate, "_")]
#      
#      output[, replicate := do.call(paste, c(.SD, sep = "_")), .SDcols = .new_reps]
#      
#      output <- output |> 
#        reshape2::acast(formula = as.formula(paste0(structure(GI_obj), collapse = " ~ ")), 
#                        value.var = "GI", 
#                        fun.aggregate = mean)
#      
#      GI_obj@blocks$collapsed <- c(GI_obj@blocks$collapsed, .collapse)
#      
#      replicates(GI_obj) <- dimnames(guideGIs(GI_obj))[[which(structure(GI_obj) == "replicate")]]
#      
#      blocks(GI_obj) <- list(map = map_replicate_layers(replicates(GI_obj)), 
#                             all = colnames(map_replicate_layers(replicates(GI_obj))), 
#                             options = colnames(map_replicate_layers(replicates(GI_obj))), 
#                             collapsed = NULL, 
#                             chosen = NULL)
#      return(GI_obj)
#    })


setMethod("add_collapsed_layers", 
          signature = "ScreenBase", 
          function(GI_obj, to_collapse = GI_obj@blocks$all) {
            
            output <- set_names(c("nothing", to_collapse)) |> 
              map(~ GI_obj@guideGIs) |> 
              map(flatten_array, dnames = GI_obj@structure, value.name = "GI") |>
              imap(~ {
                .old_reps <- GI_obj@blocks$all
                .new_reps <- setdiff(.old_reps, .y)
                .x[, c(.old_reps) := data.table::tstrsplit(replicate, "_")]
                .x[, replicate := do.call(paste, c(.SD, sep = "_")), .SDcols = .new_reps]
              }) |> 
              map(reshape2::acast, 
                  formula = as.formula(paste0(structure(GI_obj), collapse = " ~ ")), 
                  value.var = "GI", 
                  fun.aggregate = mean)
            
            
            
            GI_obj@blocks$map <- map(output, ~ dimnames(.x)[[3]]) |> 
              map(map_replicate_layers) |>
              purrr::list_flatten(name_spec = "{outer}_collapsed_{inner}_used")
            
            GI_obj@blocks$all <- names(GI_obj@blocks$map)
            GI_obj@blocks$use <- names(GI_obj@blocks$map)
            
            GI_obj@guideGIs <- output |> 
              reduce(abind::abind, along = 3)
            
            return(GI_obj)
          }
)


## rework!
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

setMethod(
  "compute_GIs", 
  signature = "ScreenBase", 
  function(GI_obj, .blocks = GI_obj@blocks$use, FDR_method = "BH") {
    
    
    stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)
    
    stopifnot("DupCorrelation values required for all blocks!" = all(.blocks %in% names(GI_obj@dupCorrelation)))
    
    for (.b in .blocks) {
      print(stringr::str_glue("{.b}..."))
      
      if (grepl("multiplex", GI_obj@screenType)) {
        
        
        geneGIs(GI_obj)[[.b]] <- set_names(GI_obj@screen_attributes$query_genes) |> 
          map(purrr::safely(\(.g) {
            .fit <- limma::lmFit(object = GI_obj@guideGIs[.g,,names(GI_obj@blocks$map[[.b]])], 
                                 block = GI_obj@blocks$map[[.b]], 
                                 correlation = GI_obj@dupCorrelation[[.b]][[.g]])
            
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
        geneGIs(GI_obj)[[.b]] <- {
          .fit <- limma::lmFit(object = GI_obj@guideGIs[,names(GI_obj@blocks$map[[.b]])], 
                               block = GI_obj@blocks$map[[.b]], 
                               correlation = GI_obj@dupCorrelation[[.b]])
          
          .efit <- limma::eBayes(.fit)
          
          output <- list(rownames(GI_obj@guideGIs), 
                         c("GI", "pval", "FDR"))
          
          output <- array(data = NA, 
                          dim = purrr::map_int(output, length), 
                          dimnames = output)
          
          output[,"GI"] <- .fit$Amean
          output[,"pval"] <- .efit$p.value[, 1]
          output[,"FDR"] <- stats::p.adjust(.efit$p.value[, 1], method = FDR_method)
          
          output
        }}
    }
    
    return(GI_obj)
  })

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
  function(GI_obj, .block) {
    
    if (missing(.block)) {
      .block <- blocks(GI_obj)$chosen
      warning(str_glue("No block set. Defaulting to using {.block}."))
    }
    
    if (grepl("multiplex", GI_obj@screenType)) {
      output <- GI_obj@geneGIs[[.block]] |> 
        flatten_array(dnames = c("query_gene", "library_gene", "variable"), 
                      value.name = "value") |>
        data.table::dcast(formula = query_gene + library_gene ~ variable, 
                          value.var = "value")
      output <- output[, .SD, .SDcols = c("query_gene", "library_gene", "GI", "pval", "FDR")]
    }
    
    if (grepl("fixed", GI_obj@screenType)) {
      output <- data.table::data.table(pair = rownames(GI_obj@geneGIs[[.block]]), 
                                       query_gene = str_split_i(rownames(GI_obj@geneGIs[[.block]]), ";", 1), 
                                       library_gene = str_split_i(rownames(GI_obj@geneGIs[[.block]]), ";", 2), 
                                       GI_obj@geneGIs[[.block]])
    }
    return(output)
  })


setMethod(
  "dupCorrelation_df", 
  signature = "ScreenBase", 
  function(GI_obj, .block) {
    
    if (missing(.block)) {
      .block <- blocks(GI_obj)$chosen
      warning(str_glue("No block set. Defaulting to using {.block}."))
    }
    
    
    if (grepl("multiplex", GI_obj@screenType)) {
      output <- data.table::data.table(
        query_gene = GI_obj@screen_attributes$query_genes, 
        data.frame(GI_obj@dupCorrelation[[.block]]))
    }
    
    if (grepl("fixed", GI_obj@screenType)) {
      output <- data.table(data.frame(GI_obj@dupCorrelation[[.block]]))
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

