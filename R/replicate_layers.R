#' @importFrom stringr str_match
#' @export replicate_layers

replicate_layers <- function(cnames) {
  layers <- list()  
  .x <- base::unlist(base::strsplit(paste0(gsub("_", ";", cnames), collapse = ";"), ";"))
  .x <- base::sapply(list(bio_rep = "b", tech_rep = "t", guides = "g"), \(.) {base::unique(grep(pattern = ., .x, value = T))})
  .x <- .x[sapply(.x, base::length) != 0]
  
  .blocks <- matrix(data = 1, nrow = base::length(cnames), ncol = base::length(.x), dimnames = list(cnames, names(.x)))
  
  .rgx <- list(bio_rep = "(b\\d+)", tech_rep = "(t\\d+)", guides = "(g\\d+_g\\d+)")
  
  for (.n in names(.rgx)) {
    if (.n %in% names(.x)) {
      .blocks[, .n] <- factor(str_match(cnames, .rgx[[.n]])[,2])
    }}
  return(.blocks)
}
