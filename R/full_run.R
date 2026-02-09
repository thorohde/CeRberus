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
  
  if (instr$verbose) {print("(1/5)")}
  
  .data <- collect_all_layer_configurations(.data)
  
  #.plt_data <- imap(.data, ~ {
  #  data.frame(config = .y, 
  #             dupcor = mean(dupCorrelation(.x), na.rm = T))}) |>
  #  rbindlist()
  
  if (!"keep_all_configurations" %in% names(instr) || (!instr$keep_all_configurations)) {
    .data <- find_optimal_configuration(.data#, verbose = instr$verbose
    )
  }
  

  #.plt_data$kept <- .plt_data$config %in% names(.data)
  
  #dupcor_plot(.plt_data, fname = file.path(instr$output_directory, "duplicateCorrelationPlot.png"))
  
  if (instr$verbose) print("(2/5)")
  
  #  if (any(map_chr(.data, is) == "PosAgnMultiplexScreen") & instr$verbose) {
  #    message("Computing symmetric GIs. This might take a while.")}
  
  .data <- map(.data, compute_models)
  

  .data <- map(.data, collect_GIs)
  
  
  if (instr$verbose) print("(4/5)")
  
  print(ls.str(.data))
  
  
  if ("output_directory" %in% names(instr) & instr$overwrite_output) {
    
    saveRDS(.data, file.path(instr$output_directory, "all_GI_objects.rds"))
    
    .output <- list()
    
    #.output$duplicate_correlation <- .plt_data
    
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

  if (instr$verbose) {print("(5/5)")}
  if (!return_output) {.data <- NULL}
  return(.data)
}