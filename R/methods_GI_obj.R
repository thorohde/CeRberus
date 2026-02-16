

.f1 <- function(x) {return(slot(x, .x))}
.f2 <- function(x, value) {slot(x, .x) <- value; return(x)}

# guideGIs methods

#purrr::walk(c("guideGIs"), ~ {
#  .x2 <- paste0(.x, "<-")
#  setMethod(.x, "guideGIs", .f1))
#  setMethod(.x2, "guideGIs", .f2)
#}
#)

# ScreenBase methods

purrr::walk(c(
  #"blocks", 
  "checks", "dupCorrelation", 
  "errors", "geneGIs", "guideGIs", "limma_models", 
  #"replicates", 
  "screen_attr", "symmGeneGIs"), ~ {
    
    .x2 <- paste0(.x, "<-")
    setMethod(.x, "ScreenBase", .f1)
    setMethod(.x2, "ScreenBase", .f2)
  }
)



setMethod("import_scores", 
          signature = signature(GI_obj = "ScreenBase"), 
          function(GI_obj) {
            
            .md <- GI_obj@metadata
            
            #input <- data.table::copy(.m$input)
            
            stopifnot("The input object needs to be a data frame." = data.table::is.data.table(.md$input), 
                      "The query gene column is not in the provided dataset." = .md$query_col %in% colnames(.md$input), 
                      "The library gene column is not in the provided dataset." = .md$lib_col %in% colnames(.md$input))
            
            setnames(.md$input, 
                     old = c(.md$bio_rep_col, .md$tech_rep_col, .md$guide_col, .md$query_col, .md$lib_col, .md$gi_col), 
                     new = c("bio_rep", "tech_rep", "guide_pair", "query_gene", "library_gene", "GI"), 
                     skip_absent = T)
            
            .md$input[, gene_pair := paste0(get("query_gene"), ";", get("library_gene"))]
            
#            .md$input[, replicate := do.call(paste, c(.SD, sep = "_")), 
#                      .SDcols = intersect(c("bio_rep", "tech_rep", "guide_pair"), colnames(.md$input))]
            
            GI_obj@metadata <- .md
            
            return(GI_obj)})





setMethod("get_screen_attributes", 
          signature = signature(GI_obj = "ScreenBase"), 
          function(GI_obj) {
            
            .i <-GI_obj@metadata$input
            .a <- list(contrasts = NULL)
            
            #pair <- NULL # to prevent package environment errors
            
            .a$query_genes <- .i[, unique(get("query_gene"))]
            .a$library_genes <- .i[, unique(get("library_gene"))]
            .a$all_genes <- union(.a$query_genes, .a$library_genes)
            .a$query_genes_not_in_lib <- setdiff(.a$query_genes, .a$library_genes)
            .a$library_genes_not_in_query <- setdiff(.a$library_genes, .a$query_genes)
            
            .a$n_query_genes <- length(.a$query_genes)
            .a$n_lib_genes <- length(.a$library_genes)
            .a$n_all_genes <- length(.a$all_genes)
            
            .a$observations_per_query <- purrr::map_int(purrr::set_names(.a$query_genes), \(.g) {.i[query_gene == .g, .N]})
            
            .a$all_pairs <- .i[, unique(get("gene_pair"))]
            .a$unique_pairs <- .i[, unique(sort_gene_pairs(get("query_gene"), get("library_gene")))]
            
            GI_obj@screen_attr <- .a
            return(GI_obj)
          })


setMethod("run_checks", 
          signature = signature(GI_obj = "ScreenBase"), 
          function(GI_obj,
                   min_query_size = 50, 
                   min_library_size = 50) {
            
            .a <- GI_obj@screen_attr
            
            GI_obj@checks <- list(
              gene_sets_equal = (length(.a$query_genes_not_in_lib) <= 0.02*.a$n_query_genes) & 
                (length(.a$lib_genes_not_in_query) <= 0.02*.a$n_lib_genes), 
              query_sufficient = .a$n_query_genes >= min_query_size, 
              library_sufficient = .a$n_lib_genes >= min_library_size, 
              stable_library_size = sum(.a$obs_per_query != median(.a$obs_per_query, na.rm = T)) <= 10, 
              sufficient_tests_per_query = sum(.a$obs_per_query >= min_library_size) >= 0.95 * length(.a$obs_per_query), 
              avg_tests_per_query = median(.a$obs_per_query, .na.rm = T)
            )
            return(GI_obj)
          })

setMethod("set_screenType", 
          signature = signature(GI_obj = "ScreenBase"), 
          function(GI_obj) {
            .md <- GI_obj@metadata
            .type <- "unknown"
            
            .c <- GI_obj@checks
            
            if (.c$gene_sets_equal & 
                .c$query_sufficient & 
                .c$library_sufficient & 
                #.attr$checks$stable_library_size & 
                .c$sufficient_tests_per_query
            ) {.type <- "MultiplexScreen"} # SymmMultiplexScreen
            
            if (!.c$gene_sets_equal & 
                .c$library_sufficient & 
                .c$query_sufficient & 
                #.attr$checks$stable_library_size & 
                .c$sufficient_tests_per_query
            ) {.type <- "MultiplexScreen"} # AsymmMultiplexScreen
            
            if (!.c$library_sufficient | 
                !.c$stable_library_size | 
                !.c$sufficient_tests_per_query
                # | #avg_tests_per_query <= 50
            ) {.type <- "FixedPairScreen"}
            
            
            # ...
            
            if (.type == "unknown") {
              warning("Unknown screen design! Forcing fixed pair run.")
              .type <- "FixedPairScreen"
            }
            
            if (.md$force_fixed_pair) {
              warning("Set up to use fixed pair structure.")
              .type <- "FixedPairScreen"
            }
            
            GI_obj <- as(object = GI_obj, Class = .type)
            
            
            if (is(GI_obj)[1] %in% c("AsymmMultiplexScreen", 
                                     "SymmMultiplexScreen", 
                                     "MultiplexScreen", 
                                     "PosAgnMultiplexScreen")) {
              GI_obj@guideGIs@space <- c("query_gene", "library_gene")
            } else if (is(GI_obj)[1] == "FixedPairScreen") {
              GI_obj@guideGIs@space <- "gene_pair"
            }
            
            GI_obj@guideGIs@replicates <- intersect(c("guide_pair", "tech_rep", "bio_rep"), 
                                                    colnames(.md$input))
            
            GI_obj@metadata <- .md
            return(GI_obj)
          })


setMethod("compute_dupCorrelation", 
          signature = signature(.x = "ScreenBase"), 
          function(.x) {
            .x@dupCorrelation <- compute_dupCorrelation(.x@guideGIs)
            return(.x)
          })











setMethod("compute_models", 
          signature = signature(GI_obj = "FixedPairScreen"), 
          function(GI_obj) {
            
            GI_obj@limma_models <- limma::lmFit(
              object = GI_obj@guideGIs@data, 
              block = GI_obj@guideGIs@blocks, 
              correlation = GI_obj@dupCorrelation) |>
              limma::eBayes()
            return(GI_obj)
          })

setMethod(
  "compute_models", 
  signature = signature(GI_obj = "MultiplexScreen"), 
  function(GI_obj) {
    output <- GI_obj@screen_attr$query_genes |>
      set_names()
    output <- output |> 
      map(safely(\(.g) {
        
        .fit <- limma::lmFit(
          object = GI_obj@guideGIs@data[.g,GI_obj@screen_attr$library_genes,], 
          block = GI_obj@guideGIs@blocks, 
          correlation = GI_obj@dupCorrelation[[.g]]) |> 
          limma::eBayes()
        
        return(.fit)
      }))
    
    GI_obj@limma_models <- map(output, "result")
    GI_obj@errors$query_genes_not_usable <- map(output, "result") |> keep(is.null) |> names()
    GI_obj@errors$GI_computation_errors <- map(output, "error")
    
    #  map("error") |> 
    #  compact() |> 
    #  map_chr(~ .x$message)
    
    if (length(GI_obj@errors$query_genes_not_usable) > 0) {
      warning(str_c("Failed computing GIs for ", length(GI_obj@errors$query_genes_not_usable), " genes."))
    }
    return(GI_obj)
  })


###











setMethod("collect_GIs", 
          signature = signature(GI_obj = "FixedPairScreen"), 
          function(GI_obj, FDR_method = "BH") {
            stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)
            
            output <- list(rownames(GI_obj@guideGIs), c("GI", "pval", "FDR"))
            
            output <- array(data = NA, 
                            dim = purrr::map_int(output, length), 
                            dimnames = output)
            
            output[,"GI"] <- GI_obj@limma_models$coefficients[, 1]
            output[,"pval"] <- GI_obj@limma_models$p.value[, 1]
            output[,"FDR"] <- stats::p.adjust(GI_obj@limma_models$p.value[, 1], method = FDR_method)
            
            geneGIs(GI_obj) <- output
            
            return(GI_obj)
          })


setMethod("collect_GIs", 
          signature = signature(GI_obj = "MultiplexScreen"), 
          function(GI_obj, FDR_method = "BH") {
            stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)
            
            output <- GI_obj@limma_models |> 
              imap(\(.m, .y) {
                .x <- data.frame(library_gene = GI_obj@screen_attr$library_genes, 
                                 query_gene = .y)
                if (is.null(.m)) {
                  .x$GI <- NA
                  .x$pval <- NA
                  .x$FDR <- NA}
                else {
                  .x$GI <- .m$coefficients[, 1]
                  .x$pval <- .m$p.value[, 1]
                  .x$FDR <- p.adjust(.m$p.value[, 1])
                }
                return(.x)})
            
            output <- output |>
              data.table::rbindlist(fill = T) |>
              data.table::melt.data.table(measure.vars = c("GI", "pval", "FDR")) |>
              reshape2::acast(formula = as.formula("query_gene ~ library_gene ~ variable"), 
                              value.var = "value", drop = F)
            
            GI_obj@geneGIs <- output
            
            if (length(setdiff(rownames(GI_obj@guideGIs), rownames(GI_obj@geneGIs))) != 0 | 
                length(setdiff(colnames(GI_obj@guideGIs), colnames(GI_obj@geneGIs))) != 0) {
              
              warning("Some genes were lost!")
              print(str(GI_obj@guideGIs))
              print(str(GI_obj@geneGIs))
              
              print(setdiff(rownames(GI_obj@guideGIs), rownames(GI_obj@geneGIs)))
              print(setdiff(colnames(GI_obj@guideGIs), colnames(GI_obj@geneGIs)))
              
            }
            return(GI_obj)
          })





setMethod(
  "collect_GIs", 
  signature = signature(GI_obj = "PosAgnMultiplexScreen"), 
  function(GI_obj, FDR_method = "BH") {
    
    GI_obj <- methods::callNextMethod(GI_obj)
    
    .x <- data.table(pair = GI_obj@screen_attr$unique_pairs)
    
    .x[, query_gene := str_split_i(pair, ";", 1)]
    .x[, library_gene := str_split_i(pair, ";", 2)]
    .x[, GI := gather_symmetric_scores(pairs = pair, .arr = GI_obj@geneGIs[,,"GI"])]
    .x[, GI_z := z_transform(GI)]
    .x[, pval := gather_symmetric_scores(pairs = pair, .arr = GI_obj@geneGIs[,,"pval"])]
    .x[, FDR := balanced_FDR(pairs = pair, 
                             pval_array = GI_obj@geneGIs[,,"pval"], 
                             fdr_method = FDR_method)]
    
    GI_obj@symmGeneGIs <- .x
    
    return(GI_obj)
  })



setMethod("screenReport", 
          signature = signature(GI_obj = "ScreenBase"), 
          function(GI_obj) {
            .a <- GI_obj@screen_attr
            .checks <- GI_obj@checks
            
            .t <- c("FixedPairScreen" = "fixed pair", 
                    "MultiplexScreen" = "multiplex", 
                    "PosAgnMultiplexScreen" = "position-agnostic multiplex")
            .t <- .t[class(GI_obj)]
            
            cat(paste0("A ", .t, " CRISPR screen with ", 
                       .a$n_query_genes, " x ", .a$n_lib_genes, " genes.\n"))
            
            #            
            #if (!.x@use_blocks) {warning("Not using any block structure.")}
            
            cat(paste0("Using ", GI_obj@guideGIs@block_layer), " as blocks. ")
            cat("")
            #print(.checks)
          })













setMethod("GI_df", 
          signature = signature(GI_obj = "FixedPairScreen"), 
          function(GI_obj) {
            output <- GI_obj@geneGIs
            
            output <- data.table::data.table(
              pair = rownames(output), 
              query_gene = str_split_i(rownames(output), ";", 1), 
              library_gene = str_split_i(rownames(output), ";", 2), 
              output)
            return(output)
          })

setMethod(
  "GI_df", 
  signature = signature(GI_obj = "MultiplexScreen"), 
  function(GI_obj) {
    
    output <- GI_obj@geneGIs
    
    output <- output |> 
      flatten_array(dnames = c("query_gene", "library_gene", "variable"), 
                    value.name = "value") |>
      data.table::dcast(formula = query_gene + library_gene ~ variable, 
                        value.var = "value")
    output[, pair := str_c(query_gene, ";", library_gene)]
    output <- output[, .SD, .SDcols = c("pair", "query_gene", "library_gene", "GI", "pval", "FDR")]
    
    return(output)
  })


setMethod(
  "GI_df", 
  signature = signature(GI_obj = "PosAgnMultiplexScreen"), 
  function(GI_obj) {
    symmGeneGIs(GI_obj)
  })



setMethod(
  "dupCorrelation_df", 
  signature = signature(GI_obj = "ScreenBase"), 
  function(GI_obj) {
    output <- data.table(data.frame(dupcor = GI_obj@dupCorrelation))
    return(output)
  })

setMethod(
  "create_log", 
  signature = signature(GI_obj = "ScreenBase"), 
  function(GI_obj) {
    list(stringr::str_glue("CRISPR screen with {GI_obj@screen_attr$n_query_genes} query genes and {GI_obj@screen_attr$n_lib_genes} library genes."), 
         stringr::str_glue("Identified screen design: {GI_obj@screenType}."), 
         "", 
         #stringr::str_glue("Possible replicate layers: {paste(GI_obj@blocks$options, collapse = ', ')}"), 
         #stringr::str_glue("Possible layer used for p-values: {GI_obj@blocks$chosen}"), 
         stringr::str_glue("Median duplicate correlation levels: "), 
         "", 
         stringr::str_glue("... and many other useful log information.")
    ) |> paste(collapse = "\n")
    
  })





setMethod("symmetry_test", 
          signature = signature(GI_obj = "MultiplexScreen"), 
          function(GI_obj, cutoff = 0.99) {
            
            .test <- map_lgl(set_names(GI_obj@guideGIs@replicates), \(.r) {
              .x <- GI_obj@guideGIs[,,.r]
              if (all(.x == t(.x)) || all(dplyr::near(.x, t(.x)))) {return(T)} else {
                return(cor(as.vector(.x), as.vector(t(.x)), use = "pairwise.complete.obs") >= cutoff)}})
            
            return(all(.test))  
            
          })
