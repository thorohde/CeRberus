#' @export replicate_layers

replicate_layers <- function(replicates) {
  layers <- list()  
  .x <- unlist(strsplit(paste0(gsub("_", ";", replicates), collapse = ";"), ";"))
  .x <- map(list(bio_rep = "b", tech_rep = "t", guides = "g"), \(.xx) {unique(grep(pattern = .xx, .x, value = T))})
  .x <- .x[map_int(.x, length) != 0]
  
  .blocks <- matrix(data = 1, nrow = length(replicates), ncol = length(.x), dimnames = list(replicates, names(.x)))
  
  .rgx <- list(bio_rep = "(b\\d+)", tech_rep = "(t\\d+)", guides = "(g\\d+)")
  
  for (.n in names(.rgx)) {
    if (.n %in% names(.x)) {
      .blocks[, .n] <- factor(str_match(replicates, .rgx[[.n]])[,2])
    }}
  return(.blocks)
}
