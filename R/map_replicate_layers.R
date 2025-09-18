#' @export

map_replicate_layers <- function(replicates) {
  layers <- list()
  .x <- unlist(strsplit(paste0(gsub("_", ";", replicates), collapse = ";"), ";"))
  .x <- purrr::map(list(bio_rep = "b", tech_rep = "t", guides = "g"), \(.xx) {unique(grep(pattern = .xx, .x, value = T))})
  
  .x <- .x |> purrr::keep(purrr::map_int(.x, length) != 0)
  
  .map <- matrix(data = 1, nrow = length(replicates), ncol = length(.x), dimnames = list(replicates, names(.x)))
  
  .rgx <- list(bio_rep = "(b\\d+)", tech_rep = "(t\\d+)", guides = "(g\\d+)")
  
  for (.n in intersect(names(.x), names(.rgx))) {
      .map[, .n] <- factor(str_match(replicates, .rgx[[.n]])[,2])
    }
  return(.map)
}
