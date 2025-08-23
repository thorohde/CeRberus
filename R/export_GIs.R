#' @export export_GIs

export_GIs <- function(GI_object, dupcor_object = NULL, directory = NULL) {
  
  print(str(GI_object))
  
  stopifnot("Output directory required!" = !is.null(directory))
  
  dir.create(directory, showWarnings = F, recursive = T)
  
  output <- list(
    GI = purrr::imap(GI_object, 
                     ~ data.table(query_gene = .y, 
                                  library_gene = rownames(.x), 
                                  .x)) |> data.table::rbindlist(), 
    dupcor = dupcor_object)
  
  paths <- c("GI_scores")
  
  if (!is.null(dupcor_object)) {
    paths <- c(paths, "duplicate_correlation")
  }
  
  paths <- purrr::map_chr(paths, ~ file.path(directory, paste0(.x, ".csv")))
  
  purrr::pwalk(list(output, paths), \(.o, .p) data.table::fwrite(x = .o, file = .p))
  

  
  return(NULL)
}