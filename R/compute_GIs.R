#' @importFrom data.table rbindlist
#' @importFrom limma lmFit eBayes
#' @importFrom purrr map map2
#' @importFrom stats p.adjust
#' @export compute_GIs

compute_GIs <- function(GI_object, dupcor, FDR_method) {
  
  if (attr(GI_object, "queried")) {
    output <- purrr::map(set_names(rownames(GI_object)), ~ GI_object[.x,,])
  }
  
  output <- purrr::map2(
    output, 
    dupcor[[attr(GI_object, "block_layer")]], 
    purrr::safely(~ {
      .fit <- limma::lmFit(.x, 
                           block = attr(GI_object, "replicate_layers")[,attr(GI_object, "block_layer")], 
                           correlation = .y)
      .efit <- limma::eBayes(.fit)
      
      data.frame(GI = .fit$Amean, 
                 pval = .efit$p.value[, 1], 
                 FDR = stats::p.adjust(.efit$p.value[, 1], method = FDR_method))
    }), .progress = T)
  
  
  output <- output |> 
    purrr::map("result") |> 
    purrr::imap(~ data.table(query_gene = .y, 
                             library_gene = rownames(.x), 
                             .x)) |> 
    data.table::rbindlist()
  
  return(output)
}

