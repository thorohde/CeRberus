
# ScreenBase methods

purrr::walk(c(
  "blocks", "block_description", "checks", #"compute_models", 
  "dupCorrelation", 
  "errors", "geneGIs", "guideGIs", "limma_models", 
  "replicates", "screen_attributes", "structure", "symmGeneGIs"), ~ {
    
    .x2 <- paste0(.x, "<-")
    setMethod(.x, "ScreenBase", function(x) return(slot(x, .x)))
    setMethod(.x2, "ScreenBase", function(x, value) {slot(x, .x) <- value; return(x)})
    #    setMethod(.x, "MultiplexScreen", function(x) return(slot(x, .x)))
    #    setMethod(.x2, "MultiplexScreen", function(x, value) {slot(x, .x) <- value; return(x)})
    
  }
)



setMethod("import_scores", signature = "ScreenBase", 
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
            
            .md$input[, pair := paste0(get("query_gene"), ";", get("library_gene"))]
            
            if (!is.null(.md$collapse_layers)) {
              #.md$input <- .md$input[, .SD, .SDcols = setdiff(colnames(.md$input), .md$collapse_layers)]
              group_cols <- setdiff(colnames(.md$input), "GI")
              
              .md$input <- .md$input[, .(GI = mean(GI, na.rm = T)), by = ..group_cols]
            }

            .md$input[, replicate := do.call(paste, c(.SD, sep = "_")), 
                      .SDcols = intersect(c("bio_rep", "tech_rep", "guide_pair"), colnames(.md$input))]
            
            GI_obj@metadata <- .md
            
            return(GI_obj)})

setMethod("get_screen_attributes", signature = "ScreenBase", 
          function(GI_obj) {
            
            .md <- GI_obj@metadata
            
            .a <- list(contrasts = NULL)
            
            #pair <- NULL # to prevent package environment errors
            
            .a$query_genes <- .md$input[, unique(get("query_gene"))]
            .a$library_genes <- .md$input[, unique(get("library_gene"))]
            .a$all_genes <- union(.a$query_genes, .a$library_genes)
            .a$query_genes_not_in_lib <- setdiff(.a$query_genes, .a$library_genes)
            .a$library_genes_not_in_query <- setdiff(.a$library_genes, .a$query_genes)
            
            .a$n_query_genes <- length(.a$query_genes)
            .a$n_lib_genes <- length(.a$library_genes)
            .a$n_all_genes <- length(.a$all_genes)
            
            .a$observations_per_query <- purrr::map_int(purrr::set_names(.a$query_genes), \(.g) {.md$input[query_gene == .g, .N]})
            
            .a$all_pairs <- .md$input[, unique(get("pair"))]
            .a$unique_pairs <- .md$input[, unique(sort_gene_pairs(get("query_gene"), get("library_gene")))]
            
            screen_attributes(GI_obj) <- .a
            return(GI_obj)
          })




setMethod("run_checks", signature = "ScreenBase", 
          function(GI_obj,
                   min_query_size = 50, 
                   min_library_size = 50) {
            
            .a <- screen_attributes(GI_obj)
            
            checks(GI_obj) <- list(
              gene_sets_equal = (length(.a$query_genes_not_in_lib) <= 0.02*.a$n_query_genes) & 
                (length(.a$lib_genes_not_in_query) <= 0.02*.a$n_lib_genes), 
              query_sufficient = .a$n_query_genes >= min_query_size, 
              library_sufficient = .a$n_lib_genes >= min_library_size, 
              stable_library_size = sum(.a$observations_per_query != stats::median(.a$observations_per_query, na.rm = T)) <= 10, 
              sufficient_tests_per_query = sum(.a$observations_per_query >= min_library_size) >= 0.95 * length(.a$observations_per_query), 
              avg_tests_per_query = stats::median(.a$observations_per_query, .na.rm = T)
            )
            return(GI_obj)
          })

#' @export

setMethod("set_screenType", signature = "ScreenBase", 
          function(GI_obj) {
            
            .type <- "unknown"
            
            .checks <- checks(GI_obj)
            
            if (.checks$gene_sets_equal & 
                .checks$query_sufficient & 
                .checks$library_sufficient & 
                #.attr$checks$stable_library_size & 
                .checks$sufficient_tests_per_query
            ) {.type <- "MultiplexScreen"} # SymmMultiplexScreen
            
            if (!.checks$gene_sets_equal & 
                .checks$library_sufficient & 
                .checks$query_sufficient & 
                #.attr$checks$stable_library_size & 
                .checks$sufficient_tests_per_query
            ) {.type <- "MultiplexScreen"} # AsymmMultiplexScreen
            
            if (!.checks$library_sufficient | 
                !.checks$stable_library_size | 
                !.checks$sufficient_tests_per_query
                # | #avg_tests_per_query <= 50
            ) {.type <- "FixedPairScreen"}
            
            
            # ...
            
            if (.type == "unknown") {
              warning("Unknown screen design! Forcing fixed pair run.")
              .type <- "FixedPairScreen"
            }
            
            if (GI_obj@metadata$force_fixed_pair) {
              warning("Set up to use fixed pair structure.")
              .type <- "FixedPairScreen"
            }
            
            GI_obj <- as(object = GI_obj, Class = .type)
            
            
            if (is(GI_obj)[1] %in% c("AsymmMultiplexScreen", 
                                     "SymmMultiplexScreen", 
                                     "MultiplexScreen", 
                                     "PosAgnMultiplexScreen")) {
              structure(GI_obj) <- c("query_gene", "library_gene", "replicate")
            } else if (is(GI_obj, "FixedPairScreen")) {
              structure(GI_obj) <- c("pair", "replicate")
            }
            
            guideGIs(GI_obj) <- reshape2::acast(
              data = GI_obj@metadata$input, 
              formula = as.formula(paste0(structure(GI_obj), collapse = " ~ ")), 
              value.var = "GI")
            
            replicates(GI_obj) <- dimnames(guideGIs(GI_obj))[[which(structure(GI_obj) == "replicate")]]
            
            return(GI_obj)
          })














setMethod("compute_dupCorrelation", signature = "FixedPairScreen", 
          function(GI_obj) {
            dupCorrelation(GI_obj) <- suppressWarnings(limma::duplicateCorrelation(
              object = GI_obj@guideGIs, 
              block = blocks(GI_obj)))$consensus.correlation
            return(GI_obj)
          })

setMethod("compute_dupCorrelation", signature = "MultiplexScreen", 
          function(GI_obj) {
            output <- GI_obj@screen_attributes$query_genes |>
              purrr::set_names() |> 
              purrr::map(\(.g) GI_obj@guideGIs[.g,,])
            
            output <- output |>
              purrr::map(\(.g) suppressWarnings(
                limma::duplicateCorrelation(object = .g, 
                                            block = blocks(GI_obj)))) |>
              purrr::map_dbl("consensus.correlation")
            dupCorrelation(GI_obj) <- output
            return(GI_obj)
          })





setMethod("compute_models", 
          signature = "FixedPairScreen", 
          function(GI_obj) {
            
            limma_models(GI_obj) <- limma::lmFit(
              object = guideGIs(GI_obj), 
              block = blocks(GI_obj), 
              correlation = dupCorrelation(GI_obj)) |>
              limma::eBayes()
            return(GI_obj)
          })

setMethod(
  "compute_models", 
  signature = "MultiplexScreen", 
  function(GI_obj) {
    output <- screen_attributes(GI_obj)$query_genes |>
      set_names()
    output <- output |> 
      map(safely(\(.g) {
        
        .fit <- limma::lmFit(
          object = guideGIs(GI_obj)[.g,screen_attributes(GI_obj)$library_genes,], 
          block = blocks(GI_obj), 
          correlation = dupCorrelation(GI_obj)[[.g]]) |> 
          limma::eBayes()
        
        return(.fit)
      }))
    
    limma_models(GI_obj) <- map(output, "result")
    errors(GI_obj)$query_genes_not_usable <- map(output, "result") |> keep(is.null) |> names()
    errors(GI_obj)$GI_computation_errors <- map(output, "error")
    
    #  map("error") |> 
    #  compact() |> 
    #  map_chr(~ .x$message)
    
    if (length(errors(GI_obj)$query_genes_not_usable) > 0) {
      warning(str_c("Failed computing GIs for ", length(errors(GI_obj)$query_genes_not_usable), " genes."))
    }
    return(GI_obj)
  })


###











setMethod("collect_GIs", 
          signature = "FixedPairScreen", function(GI_obj, FDR_method = "BH") {
            stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)
            
            output <- list(rownames(GI_obj@guideGIs), c("GI", "pval", "FDR"))
            
            output <- array(data = NA, 
                            dim = purrr::map_int(output, length), 
                            dimnames = output)
            
            output[,"GI"] <- limma_models(GI_obj)$coefficients[, 1]
            output[,"pval"] <- limma_models(GI_obj)$p.value[, 1]
            output[,"FDR"] <- stats::p.adjust(limma_models(GI_obj)$p.value[, 1], method = FDR_method)
            
            geneGIs(GI_obj) <- output
            
            return(GI_obj)
          })


setMethod("collect_GIs", 
          signature = "MultiplexScreen", 
          function(GI_obj, FDR_method = "BH") {
            stopifnot("Unknown FDR method provided." = FDR_method %in% p.adjust.methods)
            
            output <- limma_models(GI_obj) |> 
              imap(\(.m, .y) {
                .x <- data.frame(library_gene = GI_obj@screen_attributes$library_genes, 
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





setMethod(
  "collect_GIs", 
  signature = "PosAgnMultiplexScreen", 
  function(GI_obj, FDR_method = "BH") {
    
    GI_obj <- methods::callNextMethod(GI_obj)
    
    .x <- data.table(pair = screen_attributes(GI_obj)$unique_pairs)
    
    .x[, query_gene := str_split_i(pair, ";", 1)]
    .x[, library_gene := str_split_i(pair, ";", 2)]
    .x[, GI := gather_symmetric_scores(pairs = pair, .arr = geneGIs(GI_obj)[,,"GI"])]
    .x[, GI_z := z_transform(GI)]
    .x[, pval := gather_symmetric_scores(pairs = pair, .arr = geneGIs(GI_obj)[,,"pval"])]
    .x[, FDR := balanced_FDR(pairs = pair, 
                             pval_array = geneGIs(GI_obj)[,,"pval"], 
                             fdr_method = FDR_method)]
    
    symmGeneGIs(GI_obj) <- .x
    
    return(GI_obj)
  })



setMethod("screenReport", signature = "ScreenBase", 
          function(GI_obj) {
            .a <- screen_attributes(GI_obj)
            .checks <- checks(GI_obj)
            
            .t <- c("FixedPairScreen" = "fixed pair", 
                    "MultiplexScreen" = "multiplex", 
                    "PosAgnMultiplexScreen" = "position-agnostic multiplex")
            .t <- .t[class(GI_obj)]
            
            cat(paste0("A ", .t, " CRISPR screen with ", 
                       .a$n_query_genes, " x ", .a$n_lib_genes, " genes.\n"))
            
            #print(.checks)
          })













setMethod("GI_df", 
          signature = "FixedPairScreen", 
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
  signature = "MultiplexScreen", 
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
  "GI_df", signature = "PosAgnMultiplexScreen", 
  function(GI_obj) {
    symmGeneGIs(GI_obj)
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





setMethod("symmetry_test", signature = "MultiplexScreen", 
          function(GI_obj, cutoff = 0.99) {
  
  .test <- map_lgl(set_names(replicates(GI_obj)), \(.r) {
    .x <- guideGIs(GI_obj)[,,.r]
    if (all(.x == t(.x)) || all(dplyr::near(.x, t(.x)))) {return(T)} else {
      return(cor(as.vector(.x), as.vector(t(.x)), use = "pairwise.complete.obs") >= cutoff)}})
  
  return(all(.test))  
  
})
