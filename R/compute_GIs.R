#' @importFrom limma lmFit eBayes
#' @importFrom stats p.adjust
#' @export compute_GIs


compute_GIs <- function(.data, .block, .dupcor, FDR_method) {
  
  if (base::is.array(.data)) {
    .data <- purrr::map(rownames(.data), ~ .data[.x,,])
  }
  
  output <- purrr::map2(.data, .dupcor, purrr::safely(~ {
    .fit <- limma::lmFit(.x, block = .block, correlation = .y)
    .efit <- limma::eBayes(.fit)
    
    data.frame(GI = .fit$Amean, 
               pval = .efit$p.value[, 1], 
               FDR = stats::p.adjust(.efit$p.value[, 1], method = FDR_method))
  }), .progress = T)
  return(output)
}

