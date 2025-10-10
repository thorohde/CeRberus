#' @export

map_replicate_layers <- function(replicates) {
  layers <- list()
  .x <- unlist(strsplit(paste0(gsub("_", ";", replicates), collapse = ";"), ";"))
  .x <- purrr::map(list(bio_rep = "b", tech_rep = "t", guides = "g"), \(.xx) {unique(grep(pattern = .xx, .x, value = T))})
  
  .x <- .x |> purrr::keep(purrr::map_int(.x, length) != 0)
  
  output <- list(bio_rep = "(b\\d+)", 
                 tech_rep = "(t\\d+)", 
                 guides = "(g\\d+)")
  
  output <- output |> 
    purrr::keep(names(output) %in% names(.x)) |> 
    purrr::map(~ {.xx <- as.integer(factor(str_match(replicates, .x)[,2]))
    names(.xx) <- replicates
    .xx})
  
  return(output)
}
