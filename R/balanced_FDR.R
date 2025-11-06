#' @export

balanced_FDR <- function(pairs, pval_array, fdr_method) {
  
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

