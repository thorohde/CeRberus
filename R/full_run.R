#' @export


full_run <- function(yaml_fpath, return_output = T) {
  
  instr <- read_instructions(yaml_fpath)
  
  stopifnot("Output directory required!" = !is.null(instr$output_directory))
  
  dir.create(instr$output_directory, showWarnings = F, recursive = T)
  
  #for (.pkg in c("BiocManager")) {if (!require(.pkg, quietly = T)) {utils::install.packages(.pkg)}}
  #if (!require("limma", quietly = T)) {BiocManager::install("limma")}

  #if (instr$verbose) {print(str(instr, give.attr = F))}
  
  .data <- switch(tools::file_ext(instr$scores_file),
                  "csv" = data.table::fread(instr$scores_file), 
                  "rds" = readRDS(instr$scores_file))
  
  .data <- collect_all_layer_configurations(.data)
  
  #.data |> #map(~ {.x@guideGIs}) |> 
  #  str() |> print()
  
  #if (F) {
  
  #.data <- imap(.data, ~ {print(.y); compute_dupCorrelation(.x)})
  .data <- map(.data, compute_dupCorrelation)
  
  
  # store all GI objects
  if ("output_directory" %in% names(instr) & instr$overwrite_output) {
    saveRDS(.data, file.path(instr$output_directory, "all_GI_objects.rds"))
  }
    
    
  .keep_all <- "keep_all_configurations" %in% names(instr) && (isTRUE(instr$keep_all_configurations))
  print("Keep all:")
  print(.keep_all)
  .data <- find_optimal_configuration(.data, keep_all = .keep_all)
  

  .data <- .data |> 
      compute_dupcor_plot(.fpath = file.path(instr$output_directory, "duplicateCorrelationPlot.png"))
  

  .data <- map(.data, compute_models)
  
  .data <- map(.data, collect_GIs)


  if ("output_directory" %in% names(instr) & instr$overwrite_output) {
    
    .output <- list()
    
    .output$duplicate_correlation <- .data[[1]]@metadata$dupcor_data
    
      for (.n in names(.data)) {
        .output[[paste0("GI_scores_", .n)]] <- GI_df(.data[[.n]])
      }
      
      .output |>
        iwalk(~ {data.table::fwrite(
          x = .x, 
          file = file.path(instr$output_directory, paste0(.y, ".csv")))
        })
      
      #.log <- create_log(.data)
      
      #writeLines(.log, file.path(instr$output_directory, "limma.log"))
      
    }

  #}####
  
  #if (instr$verbose) {print("(5/5)")}
  if (!return_output) {.data <- NULL}
  return(.data)
}