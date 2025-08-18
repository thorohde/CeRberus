#' @importFrom data.table fwrite rbindlist
#' @importFrom utils menu
#' @importFrom purrr map imap
#' @export export_GIs

export_GIs <- function(GI_object, filepath, overwrite = F) {
  if (!base::dir.exists(base::dirname(filepath))) {
    base::dir.create(base::dirname(filepath), showWarnings = F, recursive = T)
  }
  
  .output <- GI_object |> 
    purrr::map("result") |> 
    purrr::imap(~ data.table(query_gene = .y, 
                             library_gene = base::rownames(.x), .x)) |>
    data.table::rbindlist()
  
  if (file.exists(filepath)) {
    if (overwrite) {
      data.table::fwrite(x = .output, file = filepath)
    } else {
      message("Keeping old file.")
    }
    #    overwrite <- utils::menu(c("Yes", "No"), title = base::paste0("Overwrite existing file '", filepath, "'?"))
  } else {
    data.table::fwrite(x = .output, file = filepath)
  }
  
  return(NULL)
}