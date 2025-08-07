#' @importFrom data.table fwrite rbindlist
#' @importFrom purrr map imap
#' @export export_GIs

export_GIs <- function(GI_object, filepath) {
  if (!base::dir.exists(base::dirname(filepath))) {
    base::dir.create(base::dirname(filepath))
  }
  
  .output <- GI_object |> 
    purrr::map("result") |> 
    purrr::imap(~ data.table(query_gene = .y, 
                             library_gene = base::rownames(.x), .x)) |>
    data.table::rbindlist()
  
  data.table::fwrite(x = .output, 
                     file = filepath)
  
  return(NULL)
}