#' @export

gather_symmetric_scores <- function(pairs, .arr, sep = ";") {
#  if (!isSymmetric(.arr)) {
#    warning(stringr::str_c("Input array is asymmetric! ABBA: ", 
#                           round(abba_cor(.arr), 3)))}

  genes1 <- str_split_i(pairs, sep, 1)
  genes2 <- str_split_i(pairs, sep, 2)

  stopifnot(length(setdiff(genes1, rownames(.arr))) == 0, 
            length(setdiff(genes2, colnames(.arr))) == 0)  
  
  return(purrr::map2_dbl(genes1, genes2, \(.g1, .g2) {.arr[.g1,.g2]}))}

